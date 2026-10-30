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
| Tier 3, real-world studies | amplicon profiles compared with paired shotgun profiles (gut, soil and vaginal in the paper; three more studies registered) |
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
`results/example/by_metric/`. Specificity needs the reference taxonomy (`--ref-dir`,
see [DATA.md](DATA.md)); without it the specificity columns are empty. To score your
own classifier, write its output as
`<results>/<mock>/<region>/<classifier>/<parameters>/taxonomy_final.tsv` with the
columns `Feature ID`, `Taxon` (`d__...; p__...; ...; s__...`) and `Confidence`.

Run a full cell (one mock, one region, all classifiers). The mock inputs ship in
`mocks/inputs/`; the reference databases are placed as described in [DATA.md](DATA.md):

```bash
taxonova models download all                          # Taxonova 26.10.1 or later
python scripts/01a_unpack_mock_inputs.py --mock-name mbarc26
python scripts/02_generate.py --mock-name mbarc26 --regions V4
python scripts/03_assign.py   --mock-name mbarc26 --regions V4 --n-jobs 2
python scripts/04_evaluate.py --mock-name mbarc26 --regions V4
```

`03_assign.py --dry-run` prints the commands without running them. The default
Taxonova configuration (`configs/classifiers_default.json`) uses the current release,
its model ids (`v4`, `v3v4`, `v5v7`, `16s`, `operon`) and its default thresholds.
Scoring uses the v2 rules by default ([docs/scoring.md](docs/scoring.md)).

## Reproducing the paper

The paper used Taxonova 1.0.0 with the v1.0.0 checkpoints and the v1 scoring rules.
`configs/classifiers_paper.json` pins that configuration: the checkpoint files with
their checksums, a threshold of 0.5 at every rank above species and 0.7 at species.
Install Taxonova 1.0.0 with `scripts/install_taxonova_paper.sh`
([docs/environments.md](docs/environments.md)).

```bash
python scripts/02_generate.py --mock-name mbarc26 --regions V4 --n-wildcard
python scripts/03_assign.py   --mock-name mbarc26 --regions V4 \
  --params-config classifiers_paper.json --results-dir results_paper/mbarc26
python scripts/04_evaluate.py --mock-name mbarc26 --regions V4 --rules v1 \
  --results-dir results_paper --out-dir results_paper/evaluation
```

The v1 truth treats an `N` in a source sequence as a wildcard, hence `--n-wildcard`
(also on `02b_generate_v5v7.py`). The exception is gut mock V5-V7, which the paper
scored on a corrected truth, included in `tests/fixtures/gutmock_V5-V7/`. A
regenerated cell need not match the published one: the simulation settings of each
published cell are not recorded here, and the reads depend on the ART, QIIME 2 and
DADA2 versions ([docs/scoring.md](docs/scoring.md)). On the classifier outputs and
truth behind the paper, `--rules v1` reproduces the published metric tables with a
maximum absolute difference of 0. Those outputs are not yet public; this repository
includes the gut mock V5-V7 cell, which `tests/test_regression_v1.py` checks. Results
are reported in the paper and are not repeated here.

## Adding a classifier or a mock

A classifier needs an entry in the classifier configs, a command builder in
`src/generate_commands.py` and a parser `src/process_<tool>_output.py` that writes
`taxonomy_final.tsv`. A mock needs its `<mock>_seqs.fasta` and `<mock>_taxa.txt` in
`mocks/inputs/` and an entry in `configs/mocks.json`. See
[CONTRIBUTING.md](CONTRIBUTING.md).

## Layout

| Path | Content |
|---|---|
| `scripts/` | numbered pipeline steps (`01` prepare to `06` runtime) and `loco/` |
| `src/` | command builders, output parsers, evaluator and scoring rules |
| `configs/` | paths, regions, mocks, classifier panels (`default`, `paper`, `sweep`) |
| `mocks/inputs/` | the mock community inputs behind the published cells |
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
> confidence-aware microbiome species-level amplicon classification. PLoS Comput Biol.
> 2026;22(10):e1014912. doi:[10.1371/journal.pcbi.1014912](https://doi.org/10.1371/journal.pcbi.1014912)

## License

MIT; see [LICENSE](LICENSE). Changes are listed in [CHANGELOG.md](CHANGELOG.md).
Questions and issues: <https://github.com/kardokhk/TaxoBench/issues>.
