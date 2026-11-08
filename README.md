# TaxoBench

[![tests](https://github.com/kardokhk/TaxoBench/actions/workflows/ci.yml/badge.svg)](https://github.com/kardokhk/TaxoBench/actions/workflows/ci.yml)

TaxoBench benchmarks taxonomic classifiers of rRNA amplicon sequences against a
reference-aware truth. Each query's truth is resolved against the reference the
classifier was given: where that reference has no name for the query's lineage, the
query is novel at that rank and naming it counts as over-classification. Two releases of
a classifier are compared with paired confidence intervals. It was built for
[Taxonova](https://github.com/kardokhk/taxonova); any classifier can be added.

## Install

```bash
git clone https://github.com/kardokhk/TaxoBench && cd TaxoBench
python -m venv .venv && . .venv/bin/activate
pip install -e ".[test]"
taxobench --help
```

Python 3.10 or later; scoring needs only numpy, pandas and scipy.

## How it works

A **tier** is a set of query sequences with their truth, written as tables in one data
model; a **run** is one classifier release on one tier's queries. `taxobench` converts
mock communities into a tier, runs or ingests classifier outputs, scores them and
compares releases:

```bash
taxobench import-mocks --data-dir data --ref-dir data/ref_dbs --out tiers/mocks
taxobench run --plan plan.json                               # classifiers, fail closed
taxobench score --tier tiers/mocks --runs runs/* --out scores
taxobench compare --a <run_id A> --b <run_id B> --tier tiers/mocks --runs runs/* \
  --plan compare_plan.json --out reports/A-vs-B
```

| | |
|---|---|
| Truth | per rank: known, unresolvable (a tie, never truncated), novel (no name in the reference), unscorable |
| Calls | named, group, no call, out of vocabulary (a name the run's reference lacks or names ambiguously) |
| Metrics | precision, recall, F1, Edgar accuracy, misclassification, under- and over-classification rates, coverage; precision-coverage curves and AURC; calibration (ECE, Brier); per-sample L1, Bray-Curtis, weighted UniFrac and Aitchison |
| Uncertainty | 95% cluster bootstrap CI for every rate; paired bootstrap, McNemar and Wilcoxon tests between releases, with a pre-declared confirmatory family and guardrails |
| Adapters | Taxonova, BLAST and VSEARCH (QIIME 2), QIIME 2 naive Bayes, IDTAXA, SPINGO, Kraken 2, MAPseq, and any `Feature ID, Taxon, Confidence` table |
| Mock communities | Mockrobiota-16S, MBARC-26, gut and vaginal mocks, nine MGnify biome mocks and their rRNA operon counterparts; regions V4, V3-V4, V5-V7, V1-V9 and the operon |

The commands, the data model, the run specs and the metrics are described in
[docs/taxobench.md](docs/taxobench.md); building the mock communities in
[mocks/README.md](mocks/README.md).

## Layout

| Path | Content |
|---|---|
| `src/taxobench/` | the package: data model, adapters, runner, metrics, statistics, comparison, command line |
| `mocks/` | mock community inputs, MGnify abundances and the scripts that build the mock cells |
| `tests/` | known-answer, statistics, data model, command line and adapter golden tests; mock builder tests in `tests/mocks/` |
| `docs/` | reference for the commands, data model, run specs and metrics |

## Tests

```bash
python -m pytest
```

## Citation

If you use TaxoBench, please cite the Taxonova paper:

> Kakabra K, Patel S, Temko A, Claesson MJ. Taxonova: deep learning enables fast and
> confidence-aware microbiome species-level amplicon classification. PLoS Comput Biol.
> 2026;22(10):e1014912. doi:[10.1371/journal.pcbi.1014912](https://doi.org/10.1371/journal.pcbi.1014912)

## License

MIT; see [LICENSE](LICENSE). Changes are listed in [CHANGELOG.md](CHANGELOG.md).
Questions and issues: <https://github.com/kardokhk/TaxoBench/issues>.
