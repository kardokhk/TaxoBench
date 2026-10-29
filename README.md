# Taxobench

Taxobench is a reproducible benchmarking framework for microbial amplicon
taxonomic classifiers. It was built to evaluate the [Taxonova](https://github.com/kardokhk/taxonova)
classifier against the established tools in the field, but it is written so that
any classifier and any community can be slotted in and compared on equal footing.

The whole idea is to put every classifier through the exact same pipeline,
on the exact same reads, with the exact same ground truth, so that the numbers
you get out are actually comparable. A lot of published amplicon benchmarks
compare tools that were each run with their own quirks, their own reference
database, and their own primer-trimming. Here everything shares one reference
backbone, one read-generation recipe, one evaluation routine, and one set of
metrics.

The benchmark spans several amplicon regions and, more importantly, three tiers
of community realism, from clean physical mocks all the way to paired
amplicon/shotgun real-world studies. That progression is the point: a classifier
that wins on a tidy 20-species mock is not necessarily the one you want on a real
soil sample.

## What gets compared

Eight classifiers run side by side (see `configs/classifiers_default.json`):

- **Taxonova** — the deep-learning amplicon classifier this framework was built for
- **BLASTn** (via QIIME 2)
- **VSEARCH** (via QIIME 2)
- **Kraken 2**
- **MAPseq**
- **SPINGO**
- **QIIME 2 naive Bayes** (`naive-bayes`)
- **IDTAXA**

Each classifier has its own per-method parameter block in the config (confidence
thresholds, percent-identity, k-mer size, and so on), and its own output parser
under `src/process_*_output.py` (plus the two R helpers for IDTAXA and the mmseqs
mock-prep step). Adding or swapping a classifier is mostly a matter of dropping in
a parser and a config entry — more on that below.

## Amplicon regions

All experiments are run per region. Primers and expected length windows live in
`configs/regions.json`:

| Region | Primers | Approx. length | Notes |
|--------|---------|----------------|-------|
| V4 | 515F / 806R (GTGCCAGCMGCCGCGGTAA / GGACTACNVGGGTWTCTAAT) | 200-350 bp | the classic short Illumina amplicon |
| V3-V4 | 341F / 805R (CCTACGGGNGGCWGCAG / GACTACHVGGGTATCTAATCC) | 400-500 bp | |
| V5-V7 | 799F / 1193R (AACMGGATTAGATACCCKG / ACGTCATCCCCACCTTCC) | 300-500 bp | plant/chloroplast-excluding amplicon, ~395 bp |
| V1-V9 | 27F / 1492R (AGAGTTTGATCMTGGCTCAG / TACGGYTACCTTGTTACGACTT) | 1300-1550 bp | full-length 16S |
| RRN | none | 4600-5300 bp | full rRNA operon; GROND/GTDB sequences used as-is, no primer extraction |

## The three evaluation tiers

This is the heart of Taxobench — increasingly realistic communities, all built
and scored the same way. The full mock matrix is defined in `configs/mocks.json`.

**Tier 1 — physical mock communities.** Known, lab-assembled mixtures with a
trusted species roster. Taxobench uses a concatenated Mockrobiota 16S panel
(`mockrobiota-16S`, evaluated on V4 / V3-V4 / V1-V9) and the MBARC-26 mock
(`mbarc26`, Singer 2016: 23 bacteria + 3 archaea across 10 phyla, evaluated on
V4 / V3-V4 / V5-V7 / V1-V9). These are the cleanest possible test: you know
exactly what went in.

**Tier 2 — biome-simulated mocks.** This is where realism comes in. For nine
environments — digestive system, lung, sputum, skin, soil, plants, marine,
ocean, and a mixed control — Taxobench samples Greengenes2 reference sequences
at biome-realistic genus abundances. Those abundances are not made up: they are
derived from MGnify by pulling the most abundant taxa per rank for each biome
(see `mgnify/` and the snapshot CSVs under `mgnify/results/`). The result is a
`biome_*` mock that has the species composition and the long-tail skew you
actually see in that environment, but with a known ground truth. Each biome also
has a full-operon counterpart (`gtdb_*`) for the RRN region, built from
GTDB/GROND sequences.

**Tier 3 — paired amplicon + shotgun real-world studies.** The most honest test:
real samples where there is no curated truth, so a paired shotgun-metagenomics
profile is used as a primer-bias-free reference. The amplicon ASVs get classified
by every tool, the shotgun side is profiled independently (woltka OGU tables),
and the two are compared per sample. Studies wired in (`configs/mocks.json`,
`real_world` block):

- **Zuo et al.** — pediatric ulcerative colitis gut, 16S V4 + shotgun (PRJNA759642)
- **Mantri et al.** — soil, 16S V3-V4 + shotgun (PRJNA717918 amplicon / PRJNA717813 shotgun)
- **Virtanen et al.** — vaginal, 16S V3-V4 + shotgun (PRJEB25778 / PRJEB43447)
- plus additional paired sets used during development: a Florida coral reef study
  (Cao et al., PRJNA910824) and two EMP500 sub-studies (forest soil and plant
  phyllosphere, PRJEB42019)

The expensive part of tier 3 (downloading reads, primer trimming, DADA2, shotgun
profiling) is done once per study; the cheap part (re-classifying the ASVs with an
updated model and re-comparing to the shotgun reference) can be rerun any time.

## Metrics

Evaluation is split between the mock tiers (`scripts/04_evaluate.py`,
`src/evaluate_results.py`) and the real-world tier (`scripts/05_realworld_evaluate.py`).

Per-rank classification quality (kingdom down to species):

- **precision, recall, F-measure** — species-level F is the headline number
- **detection rate / sensitivity** — fraction of true taxa recovered
- **specificity** — how well a tool refrains from predicting taxa that are not in
  the truth, measured against the per-region reference vocabulary (Greengenes2 for
  V4/V3-V4/V5-V7/V1-V9, GROND for RRN)
- **L1 / L2 abundance error** — how far off the predicted composition is

Community-level distances (especially for the real-world tier):

- **Jaccard** — taxon-set overlap
- **Bray-Curtis** — abundance-weighted dissimilarity
- **Aitchison** — Euclidean distance on CLR-transformed compositions, so it is
  compositional-aware

On top of the per-mock scores there is a misclassification vs over-classification
(MC/OC) analysis, run as a **leave-one-clade-out (LOCO)** experiment under
`scripts/loco/` and `src/loco/`. Standard mock and real-world tests cannot probe
over-classification, because every true taxon in the query is also in the
reference database. LOCO fixes that by holding entire clades out of the reference,
so you can measure how often a tool wrongly assigns a confident name to a genuinely
novel taxon (OC) versus how often it assigns the wrong name to a known one (MC),
as a function of the confidence threshold.

Everything in the evaluation is **bracket-aware**. Taxonova (and the consensus
logic generally) can emit a tie like `[species_a|species_b|species_c]` when a read
genuinely cannot be resolved past a group of equally-good labels. Both the mock
and the real-world evaluators parse those brackets and score the prediction as a
correct hit if the true taxon is inside the set, rather than treating the whole
string as one wrong label (`src/evaluate_results.py`, `scripts/04_evaluate.py`).
See "Scoring rules" below for how abstentions and groups are scored.

## Scoring rules (v1 and v2)

`scripts/04_evaluate.py --rules` selects how calls are scored. The same rules
apply to every classifier; no rule looks at which tool produced a label. The
code is in `src/scoring_v2.py` (v2) and in `scripts/04_evaluate.py` (v1).

| `--rules` | abstention (`...; g__X; s__`, `Unassigned`, `unclassified_*`) | group label `s__[A\|B]` | taxon identity in taxon-set metrics |
|---|---|---|---|
| `v2` (default) | no call: a false negative in P/R/F; left out of detection, sensitivity, specificity, Jaccard and L1/L2; reported as `unassigned_fraction` | correct if the truth is a member, at every rank where both are named | rank prefix + name |
| `v2-strict` | as `v2` | a multi-member group is no call at that rank and below | as `v2` |
| `v1` | counted as a taxon of its own, so it is a false-positive taxon in specificity, detection and L1/L2 | credited rank by rank (same outcome as `v2` when the truth has one species) | whole label string |
| `v1-patched` | as `v1`, but dropped from the taxon-set and L1/L2 metrics | as `v1` | as `v1` |

`v1` is the evaluator used for the Taxonova paper. It has two known defects that
`v2` corrects: an abstention lowers specificity, so a tool that abstains more
can look as if it calls more absent taxa; and a correct call whose upper ranks
are blank (SPINGO writes `d__; p__; ...; g__X; s__X y`) never equals the truth
string, so it counts as a false-positive taxon plus a missed one. Under `v2` a
tool that abstains a lot can reach specificity 1.0 with low recall, so report
specificity together with recall or F, not as a ranking on its own.

Under `v2` and `v2-strict` the species rows also carry `n_named_calls`,
`n_group_calls`, `group_share`, `group_size_mean`, `group_size_max`,
`group_truth_in`, `group_truth_out`, `strict_species_f`, `jaccard`,
`unassigned_fraction` and `tp` / `fp` / `fn`.

**Ground truth and N bases.** The truth builder (`src/process_mock.py`, step 02)
now treats an `N` in a source sequence as a mismatch. Up to the published runs
it was a wildcard, so a source sequence with a long run of `N` matched ASVs it
does not represent. Pass `--n-wildcard` to `02_generate.py` or
`02b_generate_v5v7.py` (or set `process_mock.N_WILDCARD = True`) to restore the
old rule. `src/rebuild_truth.py` rebuilds one truth file without QIIME 2:

```bash
python src/rebuild_truth.py \
  --observed data/mock_community/gutmock/V5-V7/paired/rep_seqs_gutmock.fasta \
  --expected data/mock_community/input_data/gutmock_seqs.fasta \
  --taxonomy data/mock_community/input_data/gutmock_taxa.txt \
  --out gutmock_V5-V7_truth.tsv            # add --n-wildcard for the v1 rule
```

With N as a mismatch the labels change in five of the published cells: gutmock
V5-V7 (125 of 138 ASVs; 124 deepen from genus to species), *Bilophila* genus to
*B. intracellularis* in biome_digestive_system and biome_sputum at V4 and V3-V4
(1 ASV each), and one mbarc26 whole-16S ASV that becomes `other`. Cells whose
truth was written by the generator rather than by sequence matching (RRN, and
whole 16S of the biome mocks, Mockrobiota-16S, gutmock and vaginalmock) are not
affected. The corrected gutmock V5-V7 truth is
`tests/fixtures/gutmock_V5-V7/data/mock_community/gutmock/V5-V7/trueish_taxonomy_gutmock.tsv`;
the copy in the original data archive was built with the wildcard rule.

**Reproducing the paper's numbers.** Use the v1 rules with the v1 truth: the
truth files from the data archive (DATA.md), with gutmock V5-V7 replaced by the
corrected file above (the paper scored that cell on the corrected truth).

```bash
python scripts/04_evaluate.py --rules v1 \
  --results-dir <results root> --truth-dir <data root holding the v1 truth> \
  --ref-dir data/ref_dbs --out-dir results/evaluation_v1
```

Run this way on the published classifier outputs, the evaluator reproduced all
26,400 values of the published `all_metrics.csv` with maximum absolute difference
0 (gutmock V5-V7 excluded there, because it was rescored on the corrected truth;
`tests/test_regression_v1.py` checks its 240 published values on a small fixture). `scripts/04b_evaluate_sweep.py`
still uses the v1 rules.

## Pipeline flow

The framework is deliberately modular: every step is a standalone CLI script under
`scripts/`, backed by importable helpers in `src/`. You can run the whole matrix
end to end, or rerun just the part you changed. Paths come from
`configs/pipeline_config.json`. There is also a more detailed walkthrough in
`docs/PIPELINE.md`.

1. **`01_prepare.py`** — prepare mock community input data. Scans the Mockrobiota
   source directory, optionally classifies expected sequences via mmseqs2, and
   writes `{seqs.fasta, taxa.txt}` pairs into the input-data directory. No-op if
   inputs already exist (unless `--force`). `01b_prepare_extra_mocks.py` adds the
   MBARC-26 mock from NCBI RefSeq records.

2. **`02_generate.py`** — generate the per-region mock datasets. For each
   (mock x region) it extracts the V-region by primer (or uses full-length GROND
   sequences directly for RRN), simulates Illumina reads with ART, runs the
   QIIME 2 demultiplex + DADA2 denoising path, and writes the "trueish" taxonomy
   that scoring is done against. `02b_generate_v5v7.py` and
   `02c_generate_v1v9_direct.py` are region-specific variants for V5-V7 and for
   curated long-read mocks.

3. **`03_assign.py`** — run all eight classifiers. For each
   (mock x region x classifier x parameter combination) it builds the command via
   `src/generate_commands.py` and runs them in parallel (joblib), then normalizes
   each tool's raw output into a common `taxonomy_final.tsv` via the
   `src/process_*_output.py` parsers.

4. **`04_evaluate.py`** — score predictions against the trueish taxonomy. Auto-
   discovers every (mock x region) that has both predictions and a truth file,
   computes all the per-rank and abundance metrics above, and writes a single
   consolidated tree under `results/evaluation/` (tidy long form plus per-metric
   wide tables). `--rules` picks the scoring rules (v2 by default; v1 for the
   published numbers), `--truth-dir`, `--ref-dir` and `--out-dir` override the
   configured locations. `04b_evaluate_sweep.py` produces precision/recall/F-vs-confidence
   sweep curves.

5. **`05_realworld_classify.py`** then **`05_realworld_evaluate.py`** — the
   real-world tier. The first runs all eight classifiers on each study's ASV FASTA
   (no truth, so just predictions). The second aggregates ASV counts to per-sample
   per-rank frequency vectors, does the same for the paired shotgun OGU table,
   joins amplicon and shotgun samples via each study's `paired_samples.csv`, and
   computes Jaccard / Bray-Curtis / Aitchison / precision / recall / F1 / L1 / L2
   per rank.

6. **`06_runtime.py`** — runtime and memory benchmark. Samples N in
   {100, 500, 1000, 2000, 5000} reads from a representative mock, runs each
   classifier a few times, and records wall time and peak RSS.

There are a few small reporting helpers alongside these (`build_tau_grid.py`,
`plot_confidence_sweep.py`, `v5v7_table1.py`).

## How the configs drive it

Everything is config-first. The JSON files under `configs/` are the single source
of truth:

- `pipeline_config.json`: directory layout (relative paths resolve against the
  repo root; `TAXOBENCH_<KEY>` environment variables, for example
  `TAXOBENCH_DATA_DIR`, override a key), conda environments, thread/job defaults
- `regions.json` — primers and length windows per amplicon region
- `mocks.json` — the full mock matrix (which mocks exist, which regions each is
  evaluated on, RRN counterparts, and the real-world study definitions)
- `classifiers_default.json` — which classifiers run and with what parameters
- `classifiers_sweep.json` — the parameter sweep variant
- `tau_grids.json` — per-classifier confidence-threshold grids for the sweep curves

**Adding a classifier:** add it to the `classifiers` list and the
`method_parameters_combinations` block in `classifiers_default.json`, point it at
a conda env in `method_conda_envs`, teach `src/generate_commands.py` how to build
its command line, and add a `src/process_<tool>_output.py` parser that maps its raw
output onto the common `taxonomy_final.tsv` schema.

**Adding a mock:** prepare its `{seqs.fasta, taxa.txt}` input (a `01*` script),
register it in `mocks.json` with the regions it should be evaluated on, and rerun
`02 -> 03 -> 04`. For a new biome, you would first refresh its abundance profile
with the `mgnify/` helpers.

## Data is not in this repo

The large inputs are intentionally left out because of size: the reference
databases (tens of GB), the generated mock reads, the real-world rep-seqs and
shotgun tables, and the trained models. None of that is tracked here. See
**[DATA.md](DATA.md)** for exactly what the pipeline expects, the on-disk layout,
where to download the reference databases, and the real-world study accessions.

## Taxonova models

The classifier this framework was built around lives at
[github.com/kardokhk/taxonova](https://github.com/kardokhk/taxonova). The
pretrained per-region Taxonova checkpoints are published on that repo's Releases
page; drop them into the reference-DB layout described in DATA.md to reproduce the
Taxonova rows.

## Tests

```bash
cd tests && python -m unittest test_scoring_rules test_truth_match test_regression_v1
```

The gutmock V5-V7 truth rebuild test runs only when the generated data tree is
present (`TAXOBENCH_DATA_DIR` or `data_dir`); it takes a few minutes.

## Documentation

- `docs/PIPELINE.md` — detailed per-step reference
- `docs/REALWORLD.md` — how to onboard a new paired study and how the shotgun comparison works
- `docs/REALWORLD_AMPLICON_QC.md` — primer-trimming / DADA2 QC audit for the real-world studies
- `docs/loco/README.md` — the leave-one-clade-out MC/OC protocol
- `docs/ENVIRONMENTS.md` — how to recreate the conda environments
- `docs/CHANGELOG.md` — the development changelog for the benchmark itself
- `mgnify/README.md` — the biome-abundance retriever

## License

MIT — see [LICENSE](LICENSE).

## Contact and acknowledgments

Questions: Kardokh Kakabra <kardokhkaka@gmail.com>.

Developed in the Claesson Lab.
