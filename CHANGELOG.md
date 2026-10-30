# Changelog

Changes to TaxoBench that affect how it runs or what it reports. The development log
kept while the benchmark was built for the paper is in
[docs/development-log.md](docs/development-log.md).

## 26.10.1 (2026-10-30)

Maintenance release alongside Taxonova 26.10.1. Scores under `--rules v1` and `--rules v2`
are unchanged.

**Taxonova configuration**

- `configs/classifiers_default.json` runs the current Taxonova release
  (`taxonova>=26.10.1`) through `taxonova classify`, with the model ids `v4`, `v3v4`,
  `v5v7`, `16s` and `operon` and the package's default thresholds (species 0.8).
- `configs/classifiers_paper.json` pins the published runs: Taxonova 1.0.0, the v1.0.0
  checkpoints (path, URL, md5 and sha256), 0.5 at every rank above species, 0.7 at
  species. Score with `--rules v1`. `classifiers_sweep.json` uses the same settings.
- A pinned checkpoint is checked against its md5 before use.
- Behavior change for custom configs without a `taxonova` block: the parent-rank
  thresholds are no longer forced to 0.5; the package defaults apply. Add
  `"parent_gates"` to keep the old behavior.

**Fixes ported from the benchmark used for the paper**

- QIIME 2 query artifacts are written per classifier, so concurrent BLAST, VSEARCH and
  naive Bayes runs no longer overwrite each other; an artifact older than its FASTA is
  imported again.
- `05_realworld_classify.py` checks that every output covers exactly the query set.
- LOCO: Kraken 2 and MAPseq in the panel; replicate roots and tags
  (`--loco-root`, `--tag`, `--fixed-higher-ranks`); a non-empty hold-out is never
  overwritten; MAPseq `N/A`, `no blast hit` and `ambiguous` count as no call;
  `04c_evaluate_opengates.py`; VSEARCH `--maxrejects`.
- `build_tau_grid.py`: MAPseq treated as continuous; dense percentile grid.

**Other fixes**

- Scripts no longer depend on `/path/to` placeholders or cluster paths; all roots come
  from `configs/pipeline_config.json`. `02b_generate_v5v7.py` failed to import on any
  checkout.
- `scripts/02d_generate_longread_biomes.py`: the long-read biome generator, previously
  referenced in a source tree outside this repository. Its seed no longer depends on
  Python's salted string hash.
- `mgnify/results/2025-08-20_snapshot/`, the MGnify tables behind the biome mocks, is
  now tracked.
- `configs/mocks.json` lists the gut and vaginal mocks and the V5-V7 region.
- Model globs are sorted; `--dry-run` creates no directories; `generate_commands`
  imports QIIME 2 only when needed; the Taxonova training fallback had an invalid
  output path and flag; `loco/03_classify.py --help` failed.
- Environment paths are set once, in `conda_envs` of `pipeline_config.json`;
  `conda_bin` or `TAXOBENCH_CONDA_BIN` sets the conda executable.

**Repository**

- Shorter README; design, scoring, pipeline and environment notes in `docs/`
  (lower-case file names, as cited in the paper).
- Continuous integration (`.github/workflows/ci.yml`), `requirements-eval.txt`,
  `CITATION.cff`, `CONTRIBUTING.md`, tests for the configurations and command building.

## 2026-10-29: scoring v2

- `scripts/04_evaluate.py --rules {v2,v2-strict,v1,v1-patched}`, `v2` by default:
  an abstention is no call; a group label is correct if the truth is a member; taxon
  identity is rank and name. `v1` reproduces the published values. New options
  `--truth-dir`, `--ref-dir`, `--out-dir`.
- Truth builder: N in a source sequence is a mismatch (`--n-wildcard` restores the
  published rule); `src/rebuild_truth.py`.
- Relative paths in `pipeline_config.json`; `TAXOBENCH_<KEY>` overrides.
- Tests, including a regression against the 240 published gut mock V5-V7 values.

Details: [docs/scoring.md](docs/scoring.md).
