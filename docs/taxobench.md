# TaxoBench reference

The `taxobench` package answers one question per release: is classifier release B
better than release A, on which tier, region, reference and rank, by how much and with
what uncertainty.

```bash
pip install -e ".[workflow,test]"
taxobench --help
taxobench adapters                       # installed adapters and what their confidence means
```

## Commands

| command | does |
|---|---|
| `import-mocks` | converts the mock builders' output (`mocks/`) into a tier with query FASTA files |
| `run` | runs classifiers from one spec or a plan, fail closed, and writes their predictions |
| `ingest` | parses an existing raw output into prediction tables, with the same checks |
| `validate` | checks tier and run tables against the data model |
| `score` | scores runs against a tier: metrics M1 to M11 with 95% bootstrap CIs |
| `compare` | paired comparison of two runs, with a confirmatory family and guardrails |
| `report` | Markdown summary of a `scores.tsv` |
| `adapters` | lists the parsers and their confidence kinds |

Every command exits with a non-zero status on any failure.

## Data model

Seven long tables, tab-separated with a header, each file with a JSON sidecar holding
its sha256 (`<file>.json`, or a `MANIFEST.json` in the directory). A table may be split
into `<table>.<part>.tsv` files. Extra columns are kept and ignored. A tier is a
directory with the first five tables; a run directory holds the last two.

| table | required columns | key |
|---|---|---|
| `references` | ref_id, release | ref_id |
| `taxa` | ref_id, rank, taxon_id, name, parent_taxon_id (optional is_group) | ref_id, taxon_id |
| `queries` | query_id, tier, dataset, region, split (`calib`/`eval`), cluster_id (or `cluster`) | query_id |
| `abundance` | query_id, sample_id, count | query_id, sample_id |
| `truth` | query_id, ref_id, rank, truth_set, status | query_id, ref_id, rank |
| `runs` | run_id, classifier, version, ref_id | run_id |
| `predictions` | run_id, query_id, rank, raw_label, taxon_ids, confidence, call | run_id, query_id, rank |

- Ranks: domain, phylum, class, order, family, genus, species. Truth has a row for every
  rank of every query; predictions likewise.
- Truth status: `known` (one taxon), `unresolvable` (a tie of several taxa, never
  truncated), `novel` (the reference has no name for this lineage at the rank),
  `unscorable` (excluded).
- Prediction call: `named`, `group`, `nocall`, `oov` (a label the run's reference does
  not contain, or names ambiguously).
- `truth_set` and `taxon_ids` join taxon_ids with `|`. A taxon flagged `is_group` stands
  for its bracketed members when they exist as taxa of the same rank.

Any tier written in this format can be scored: held-out genomes, real-read mocks, or
the mock communities of `mocks/`.

## Mock communities

The builders in `mocks/scripts/` write each cell's truth table and query sequences under
`data/mock_community/<mock>/<region>/` ([mocks/README.md](../mocks/README.md)).
`import-mocks` turns them into a tier:

```bash
taxobench import-mocks --data-dir data --ref-dir data/ref_dbs --out tiers/mocks
taxobench import-mocks --data-dir data --ref-dir data/ref_dbs --out tiers/mbarc26 \
  --cell mbarc26:V4 --cell mbarc26:V1-V9
```

It writes the tier tables and one FASTA per cell in `tiers/mocks/queries/`, with query
IDs `<mock>__<region>__<id>`, so classifier outputs on those files score directly.
References are Greengenes2 2024.09 (`<ref-dir>/<region>/taxonomy.tsv`, reference ids
`gg2-2024.09/<region>`) and GROND release 207 for the operon
(`<ref-dir>/RRN/grond_clean.tsv`, `grond-r207/RRN`). A truth name the region's reference
does not contain is `novel` at that rank. Every query is in split `eval` and is its own
cluster, and the abundance table is empty, so the per-sample profile metrics (M10, M11)
are not computed for mock tiers.

## Adapters

One per tool, versioned, with golden tests on raw outputs of the real binaries
(`tests/test_adapters_golden.py`). Every label is resolved against the run's
reference: a name resolves when exactly one taxon of that rank has it (after narrowing by
the resolved parent). A repair proposed by an adapter is accepted only when it, too,
leads to one taxon:

- SPINGO writes `Enterococcus J nangangensis` for `Enterococcus_J nangangensis`; the
  reference name equal after `_` -> space is used when unique. Ranks above the genus come
  from the reference.
- Kraken 2 databases can drop the numeric suffix of a GTDB-style genus; reference names
  with suffixes removed are matched when unique, or the genus is settled by a fully
  named species below it. Its per-taxon value (report percentage / 100) is the clade's
  share of reads in the run, kept as `run_share`: not a per-read value, not a probability.
- MAPseq: the first rank below the cutoff ends the lineage.

| adapter | input | confidence |
|---|---|---|
| taxonova | `taxonova classify` (or `classification`) table | probability |
| idtaxa | output of `taxobench/tools/run_idtaxa_classification.R` | score |
| spingo | SPINGO 1.3 output | score |
| kraken2 | `--output` with `--use-names`, plus `--report` | run_share |
| mapseq | `.mseq` | score |
| blast, vsearch | q2-feature-classifier table | score |
| naive-bayes | q2-feature-classifier table | probability |
| taxonomy-tsv, spingo-final, kraken2-final | a `Feature ID, Taxon, Confidence` table (precomputed outputs) | as named |

Third-party adapters register under the `taxobench.adapters` entry point.

## Running classifiers

A run spec names the adapter, version, reference, query FASTA, command (argv, no shell)
and raw output. `taxobench run --spec` clears the output directory, runs the command
(locally, in an Apptainer image or a conda prefix), checks exit status, freshness and
query coverage, and only then writes `predictions.tsv`, `runs.tsv` and `run.json`.

```json
{
  "adapter": "taxonova", "classifier": "taxonova", "version": "26.11.2",
  "ref_id": "gg2-2024.09/V4", "taxa": "tiers/mocks",
  "query_fasta": "tiers/mocks/queries/mbarc26__V4.fasta",
  "command": ["taxonova", "classify", "{query}", "--model", "v4", "-o", "{outdir}/raw.tsv"],
  "raw": "raw.tsv",
  "runner": {"type": "conda", "prefix": "/path/to/env"}
}
```

`{query}`, `{outdir}`, `{threads}`, `{model}` and the keys of `parameters` are filled in
before the command runs. `taxobench ingest` does the same checks for an existing raw
output. `taxobench run --plan plan.json` runs a plan (`{"defaults": {...}, "runs":
[spec, ...]}`) one spec after another; `--snakemake --cores N` schedules the specs in
parallel with Snakemake (`pip install -e ".[workflow]"`). On a cluster, running
Snakemake inside one node allocation packs many small runs into its cores.

## Metrics

Per run, cell (tier x dataset x region x reference) and rank, on `eval` queries.
K = known or unresolvable rows, L = novel rows, unscorable rows excluded.

- M1 outcomes: TP, MC (named, inconsistent at this or a higher rank), UC (not named;
  `uc_wrong_above` counts those wrong higher up), OC (named on a novel row), CA.
- M2 to M6: precision TP/(TP+MC+OC) (NaN when nothing is named), recall TP/K, F1 from
  unrounded counts, Edgar accuracy TP/(K+OC), MCR, UCR, OCR, coverage.
- M7 group rules: G1 group-containing (default), G2 strict (always reported beside it),
  G3 size-weighted.
- M8 precision-coverage curve, step-wise average precision, AURC, F1-optimal threshold
  fitted on `calib` and applied to `eval` (probability and score confidences only).
- M9 calibration: ECE with 15 equal-mass bins (ties kept together) and 10 equal-width
  bins, Brier score (probability confidences only).
- M10 per-sample read-weighted profiles: L1 and Bray-Curtis under convention A
  (renormalised over named mass) and B (unassigned kept); weighted UniFrac on the
  reference taxonomy with unit branch lengths. Group mass sits at the members' lowest
  common ancestor and is never resolved with the truth.
- M11 Aitchison over truth taxa, one pooled "other named" part and "unassigned", with a
  pseudocount of 0.5 read.

Every rate has a 95% cluster bootstrap CI (`--n-boot`, seed recorded).

## Comparing releases

```bash
taxobench compare --a <run_id A> --b <run_id B> --tier tiers/T1 tiers/T2 --runs runs/ \
  --plan plan.json --out reports/A-vs-B
```

Both runs must use the same reference and have predicted the same queries in every cell.
Per cell, rank, rule and metric: A, B, delta = B - A with a paired cluster bootstrap 95%
CI (10,000 resamples, the same clusters for both), bootstrap *P* and an exact McNemar
test. Per tier: the cells pooled, plus the median delta over cells, Wilcoxon *P* and a
Hodges-Lehmann CI. The plan declares, before the run, the confirmatory family (Holm) and
the non-inferiority guardrails; everything else is Benjamini-Hochberg adjusted and
labelled exploratory. Output: `compare.md`, `compare.tsv`, `guardrails.tsv`. The command
exits 3 when a guardrail fails.

```json
{"confirmatory": [{"tier": "T1", "rank": "species", "metric": "f1"},
                  {"tier": "T2", "rank": "species", "metric": "f1"}],
 "guardrails":   [{"tier": "T1", "rank": "genus", "metric": "f1", "margin": 0.01}]}
```
