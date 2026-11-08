# Changelog

Changes to TaxoBench that affect how it runs or what it reports. The development log
kept while the benchmark was built for the paper is in the git history
(`docs/CHANGELOG.md` at commit c9c4f46).

## 2.0.0 (2026-11-08)

TaxoBench 2: the benchmark is a Python package with one command, `taxobench`. Truth,
name resolution and metrics are defined anew, so its scores are not directly comparable
with those of 26.11.2 and earlier.

**Package**

- `pip install -e .` installs `taxobench` with the commands `import-mocks`, `run`,
  `ingest`, `validate`, `score`, `compare`, `report` and `adapters`
  ([docs/taxobench.md](docs/taxobench.md)).
- One data model of seven tables with sha256 sidecars for tiers (references, taxa,
  queries, abundance, truth) and runs (runs, predictions). Truth has a status per rank:
  known, unresolvable, novel or unscorable.
- Adapters parse each tool's raw output and resolve every name against the run's
  reference; a name is repaired only when the repair leads to exactly one taxon, and a
  name the reference lacks is `oov`. Golden tests on raw outputs of all eight tools.
- Metrics M1 to M11 under three group rules, with 95% cluster bootstrap CIs.
- `compare`: paired bootstrap CIs, exact McNemar, Wilcoxon and Hodges-Lehmann per tier,
  a pre-declared confirmatory family (Holm) and non-inferiority guardrails (exit 3 when
  one fails).
- `run`: fail-closed runner for local, conda or Apptainer commands, and Snakemake
  scheduling of a plan.
- `import-mocks`: the mock builders' output becomes a tier, with per-cell query FASTA
  files whose IDs are the tier's query IDs.

**Repository**

- The mock builders and their configs moved to `mocks/scripts/` and `mocks/configs/`.
- Removed, replaced by the package: the classifier runner, parsers and scorer scripts
  (`scripts/03_assign.py` to `06_runtime.py`, `src/`), their configs, the real-world
  and LOCO scripts, and the documentation of those scripts. The code of 26.11.2 stays
  available at tag `v26.11.2`.
- Python 3.10 to 3.12 in CI.

## 26.11.2 (2026-11-08)

Release alongside Taxonova 26.11.2. Scores of valid outputs under `--rules v1` and
`--rules v2` are unchanged: rescoring the paper's classifier outputs with `--rules v1`
matches all 47 published cells (15,040 values) with a difference of 0.

**Taxonova configuration**

- `configs/classifiers_default.json` runs `taxonova>=26.11.2` without a threshold, so
  each model applies its own calibrated per-rank cut-offs (the package default since
  Taxonova 26.11.0). Its results are written under `taxonova/default/` instead of
  `taxonova/confidence_threshold-0.8/`. `classifiers_paper.json` and
  `classifiers_sweep.json` are unchanged.

**Repository**

- Fewer files at the top level: `DATA.md` is now `docs/data.md`, `CONTRIBUTING.md` is in
  `.github/`, `requirements-eval.txt` is `requirements.txt`, and `mgnify/` is in
  `mocks/mgnify/`.
- The download-based check of all published cells and its weekly CI job are removed,
  with `requirements-golden.txt`; `tests/test_regression_v1.py` keeps checking the
  published gut mock V5-V7 cell.

**Fail closed**

- `03_assign.py` and `05_realworld_classify.py` delete the outputs of any earlier run
  in a job's directory before the job runs, so a failed rerun can no longer leave the
  previous `taxonomy_final.tsv` to be scored. A job counts as failed unless it exits 0
  and writes a new `taxonomy_final.tsv` holding every query ID exactly once; a failed
  job's output is renamed `taxonomy_final.tsv.failed`. Each job writes
  `run.status.json`. Missing references, command build errors and failed jobs make the
  script exit 1.
- `04_evaluate.py` refuses an output with missing, unexpected or duplicated IDs, one
  whose `run.status.json` reports a failure or a later start than the file, and
  predictions without a truth file. It then exits 1 and writes nothing. Finding
  nothing to score also exits 1.
- `02_generate.py` exits 1 when the QIIME 2/DADA2 step or the truth build fails; an
  `art_illumina` failure raises an error. The IDTAXA, MAPseq and Taxonova parsers exit 1
  instead of writing an empty table.
- `conda_bin` and `TAXOBENCH_CONDA_BIN` are honoured by every command builder (they were
  used only by the LOCO scripts); `micromamba` works.

**Paper reproduction**

- `scripts/install_taxonova_paper.sh` installs the dependencies of Taxonova 1.0.0,
  pinned to the paper environment (torch 2.7.1, numpy 1.26.4, pandas 2.2.2,
  scikit-learn 1.4.2, ...), then the package with `--no-deps`, and checks that
  `taxonova classification` starts. The v1.0.0 `pyproject.toml` declares no
  dependencies, so the previous script installed Taxonova without any and the command
  failed on import.

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
- `scripts/install_taxonova_paper.sh` installs Taxonova 1.0.0; the `v1.0.0` tag cannot
  be installed with `pip install git+...` because its `setup.py` names a missing
  launcher.
- With a model id, Taxonova runs even when `data/ref_dbs/<region>/` does not exist.
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

**Mock inputs**

- `mocks/inputs/` holds the 44 mock input files behind the published cells (gzipped,
  with md5 sums), and `scripts/01a_unpack_mock_inputs.py` places them. Nothing in the
  repository could rebuild the biome inputs, and `01b_prepare_extra_mocks.py` writes
  MBARC-26 labels that differ from the published ones.

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
  (lower-case file names, as cited in the paper). `docs/realworld.md` describes the
  preprocessing steps by tool; the per-study scripts are not yet public. The
  reproduction notes state that the v1 truth needs `--n-wildcard`.
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

Details: [docs/scoring.md](https://github.com/kardokhk/TaxoBench/blob/v26.11.2/docs/scoring.md) at tag `v26.11.2`.
