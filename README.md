# TaxoBench

[![tests](https://github.com/kardokhk/TaxoBench/actions/workflows/ci.yml/badge.svg)](https://github.com/kardokhk/TaxoBench/actions/workflows/ci.yml)

TaxoBench is a benchmark for taxonomic classifiers of rRNA amplicon sequences. It
runs every classifier on the same simulated and physical mock communities, with the
same reference data and the same ground truth, and scores the outputs with one
evaluator. It was built for the [Taxonova](https://github.com/kardokhk/taxonova)
paper and is maintained alongside Taxonova; any classifier can be added.

## What it provides

| | |
|---|---|
| Regions | 16S V4, V3-V4, V5-V7, whole 16S (V1-V9), rRNA operon (RRN) |
| Tier 1, physical mocks | Mockrobiota-16S, MBARC-26, gut mock, vaginal mock |
| Tier 2, biome-simulated mocks | nine biomes (digestive system, lung, sputum, skin, soil, plants, marine, ocean, mixed) with MGnify abundances; RRN counterparts from GROND |
| Tier 3, real-world studies | amplicon profiles compared with paired shotgun profiles (gut, soil, vaginal) |
| Over-classification test | leave-one-clade-out (LOCO) reference hold-outs |
| Classifiers | Taxonova, BLASTn and VSEARCH (QIIME 2), QIIME 2 naive Bayes, IDTAXA, SPINGO, Kraken 2, MAPseq |
| Metrics | precision, recall, F-measure, detection rate, specificity, L1 and L2 abundance error per rank; Jaccard, Bray-Curtis and Aitchison for the real-world tier |

Details: [docs/benchmark.md](docs/benchmark.md).

## Quick start

Scoring needs only Python 3.11 and a few packages. Running classifiers needs the
conda environments in [docs/environments.md](docs/environments.md).

```bash
git clone https://github.com/kardokhk/TaxoBench && cd TaxoBench
python -m venv .venv && . .venv/bin/activate
pip install -r requirements-eval.txt
```

Score classifier outputs against a mock. The repository includes one published cell
(gut mock, V5-V7, eight classifiers) as an example:

```bash
python scripts/04_evaluate.py \
  --results-dir tests/fixtures/gutmock_V5-V7/results \
  --truth-dir   tests/fixtures/gutmock_V5-V7/data \
  --out-dir     results/example
```

The scores are in `results/example/all_metrics.csv` (one row per mock, region,
classifier, parameter set, rank and metric) and in one table per metric under
`results/example/by_metric/`. Specificity needs the reference taxonomy
(`--ref-dir`, see DATA.md) and is reported as NaN without it. To score your own classifier, write its output as
`<results>/<mock>/<region>/<classifier>/<parameters>/taxonomy_final.tsv` with the
columns `Feature ID`, `Taxon` (`d__...; p__...; ...; s__...`) and `Confidence`.

Run a full cell (one mock, one region, all classifiers), after placing the data as
described in [DATA.md](DATA.md):

```bash
taxonova models download all                          # Taxonova 26.10.1 or later
python scripts/02_generate.py --mock-name biome_soil --regions V4
python scripts/03_assign.py   --mock-name biome_soil --regions V4 --n-jobs 2
python scripts/04_evaluate.py --mock-name biome_soil --regions V4
```

`03_assign.py --dry-run` prints the commands without running them. The default
Taxonova configuration (`configs/classifiers_default.json`) uses the current release,
its model ids (`v4`, `v3v4`, `v5v7`, `16s`, `operon`) and its default thresholds.
Scoring uses the v2 rules by default ([docs/scoring.md](docs/scoring.md)).

## Reproducing the paper

The paper used Taxonova 1.0.0 with the v1.0.0 checkpoints and the v1 scoring rules.
`configs/classifiers_paper.json` pins that configuration: the checkpoint files with
their md5 sums, a threshold of 0.5 at every rank above species and 0.7 at species.

```bash
# Taxonova as published, in its own environment (docs/environments.md)
python scripts/03_assign.py --mock-name biome_soil --regions V4 \
  --params-config classifiers_paper.json --results-dir results_paper/biome_soil

# Score with the published rules and the published truth
python scripts/04_evaluate.py --rules v1 \
  --results-dir results_paper --truth-dir <v1 truth root> --out-dir results/evaluation_paper
```

The v1 truth is the generated truth of the benchmark data, with the gut mock V5-V7
cell replaced by the corrected file in
`tests/fixtures/gutmock_V5-V7/data/mock_community/gutmock/V5-V7/`. On the classifier
outputs behind the paper, `--rules v1` reproduces the published metric tables with a
maximum absolute difference of 0. `tests/test_regression_v1.py` checks this on
the gut mock V5-V7 cell. Results are reported in the paper; they are not repeated here.

## Adding a classifier or a mock

A classifier needs an entry in the classifier configs, a command builder in
`src/generate_commands.py` and a parser `src/process_<tool>_output.py` that writes
`taxonomy_final.tsv`. A mock needs its `<mock>_seqs.fasta` and `<mock>_taxa.txt` in
the input directory and an entry in `configs/mocks.json`. See
[CONTRIBUTING.md](CONTRIBUTING.md).

## Layout

| Path | Content |
|---|---|
| `scripts/` | numbered pipeline steps (`01` prepare to `06` runtime) and `loco/` |
| `src/` | command builders, output parsers, evaluator and scoring rules |
| `configs/` | paths, regions, mocks, classifier panels (`default`, `paper`, `sweep`) |
| `mgnify/` | MGnify abundance retriever and the snapshot used for the biome mocks |
| `tests/` | unit and regression tests, with a small fixture |
| `docs/` | benchmark design, scoring rules, pipeline, environments, real-world tier, LOCO |
| `DATA.md` | reference databases, data layout and study accessions |

## Tests

```bash
cd tests && python -m unittest discover -s . -p "test_*.py"
```

One test rebuilds the gut mock V5-V7 truth and runs only when the benchmark data are
present (`TAXOBENCH_DATA_DIR`).

## Citation

If you use TaxoBench, please cite the Taxonova paper:

> Kakabra K, Patel S, Temko A, Claesson MJ. Taxonova: deep learning enables fast and
> confidence-aware microbiome species-level amplicon classification. PLOS Comput Biol.
> 2026. doi:[10.1371/journal.pcbi.1014912](https://doi.org/10.1371/journal.pcbi.1014912)

## License

MIT; see [LICENSE](LICENSE). Changes are listed in [CHANGELOG.md](CHANGELOG.md).
Questions and issues: <https://github.com/kardokhk/TaxoBench/issues>.
