# Pipeline reference

Every step is a standalone script under `scripts/`, with helpers in `src/`. Run the
whole matrix or only the part that changed. Paths come from
`configs/pipeline_config.json`; relative paths resolve against the repository root,
and `TAXOBENCH_<KEY>` environment variables (for example `TAXOBENCH_DATA_DIR`)
override a key. Every script answers `--help`.

## Steps

```
01_prepare.py, 01b_prepare_extra_mocks.py
        -> data/mock_community/input_data/<mock>_seqs.fasta, <mock>_taxa.txt
02_generate.py, 02b_generate_v5v7.py          (V4, V3-V4, V5-V7: ART and DADA2)
02c_generate_v1v9_direct.py                   (V1-V9 of curated mocks)
02d_generate_longread_biomes.py               (V1-V9 and RRN of biome mocks)
        -> data/mock_community/<mock>/<region>/.../rep_seqs_<mock>.{fasta,qza}
           data/mock_community/<mock>/<region>/trueish taxonomy
03_assign.py
        -> results/<mock>/<region>/<classifier>/<param-id>/taxonomy_final.tsv
04_evaluate.py
        -> results/evaluation/all_metrics.csv, by_metric/<metric>.csv
04b_evaluate_sweep.py, plot_confidence_sweep.py   (confidence sweeps)

05_realworld_classify.py, 05_realworld_evaluate.py   (tier 3)
06_runtime.py                                        (time and memory)
loco/01 to loco/05                                   (docs/loco/README.md)
```

`01_prepare.py` builds Mockrobiota inputs from a Mockrobiota checkout and does nothing
if the inputs exist (unless `--force`). `01b_prepare_extra_mocks.py` fetches MBARC-26
from NCBI RefSeq.

## Step 02: mocks and truth

`02_generate.py` extracts the region by primer match, simulates Illumina reads with ART,
demultiplexes and denoises with QIIME 2 and DADA2, and writes the truth by matching
each ASV to its source sequences. It refuses V1-V9 and RRN, because ART reads are too
short. For those regions, `02c_generate_v1v9_direct.py` (curated mocks) and
`02d_generate_longread_biomes.py` (biome mocks) mutate full-length sequences directly
at a per-sequence rate drawn from 0.002 to 0.02, so the truth is known by
construction. `--n-wildcard` restores the published truth rule for N bases
([scoring.md](scoring.md)).

## Step 03: classifiers

`03_assign.py` builds one command per mock, region, classifier and parameter set
(`src/generate_commands.py`) and runs them in parallel with joblib (`--n-jobs`). Each
tool's output is converted to `taxonomy_final.tsv` (`Feature ID`, `Taxon`,
`Confidence`) by `src/process_<tool>_output.py` or `src/qiime_helper.py`.
`--dry-run` prints the commands; `--params-config` selects the classifier config;
`--ref-root` selects another reference root.

| Classifier | Command | Reference input under `data/ref_dbs/<region>/` |
|---|---|---|
| taxonova | `taxonova classify` (default) or `taxonova classification` (paper) | model id, or a checkpoint (`taxonova/*model.pt`) |
| vsearch | `qiime feature-classifier classify-consensus-vsearch` | `vsearch/{sequences,taxonomy}.qza` |
| blast | `qiime feature-classifier classify-consensus-blast` | `blast/{blastdb,taxonomy}.qza` |
| naive-bayes | `qiime feature-classifier classify-sklearn` | `naive-bayes/*sklearn*.qza` |
| idtaxa | `Rscript src/run_idtaxa_classification.R` | `idtaxa/*.rds` |
| spingo | `spingo` (own environment) | `spingo/*.fasta` |
| kraken2 | `kraken2` | `kraken2/{hash,opts,taxo}.k2d` |
| mapseq | `mapseq` | `mapseq/reference.{fasta,tax}` |

### Classifier configs

| File | Purpose |
|---|---|
| `classifiers_default.json` | current Taxonova release (`taxonova>=26.10.1`), model ids, package default thresholds |
| `classifiers_paper.json` | the published runs: Taxonova 1.0.0, v1.0.0 checkpoints pinned by md5, 0.5 above species, 0.7 at species |
| `classifiers_sweep.json` | confidence sweep for `04b_evaluate_sweep.py`, with the paper's Taxonova settings |

Each config lists the classifiers, their environments (`method_conda_envs`, keys of
`conda_envs` in `pipeline_config.json`) and their parameter grids
(`method_parameters_combinations`; every combination is one run, named by its
parameters). The `taxonova` block sets:

| Key | Meaning |
|---|---|
| `requirement` | the Taxonova package the environment must provide |
| `install` | how to install it, where `pip install` of the requirement is not enough |
| `cli` | `classify` (model ids; taxonova 26.10.1 or later) or `classification` (every release) |
| `models` | per region: a model id (`v4`), a checkpoint path relative to the reference root, or `{path, url, md5, sha256}` |
| `parent_gates` | `null` for the package defaults, or a threshold per rank above species |

Without a `models` entry the first `taxonova/*model.pt` of the region is used.

Default parameters:

| Classifier | Parameter | Default panel |
|---|---|---|
| taxonova | `confidence_threshold` (species) | 0.8 (paper: 0.7) |
| vsearch, blast | `p-perc-identity`, `p-min-consensus`, `p-maxaccepts` | 0.80, 0.51, 10 |
| naive-bayes | `p-confidence` | 0.7 |
| idtaxa | `p-confidence` (0 to 100) | 60 |
| spingo | `kmersize`, `bootstrap`, `subsample` | 8, 10, 10 |
| kraken2 | `confidence` | 0.0 |
| mapseq | `cutoff` | 0.4 |

## Step 04: evaluation

`04_evaluate.py` finds every `<mock>/<region>` that has predictions and a truth file
and writes:

- `all_metrics.csv`: long format, columns `mock, region, method, param, level, metric,
  value`;
- `by_metric/<metric>.csv`: one table per metric, rows `mock, region, method, param`,
  columns `Class, Order, Family, Genus, Species`.

The output directory is replaced as a whole on each run; predictions are never
modified. Options: `--rules {v2,v2-strict,v1,v1-patched}` (default `v2`),
`--truth-dir`, `--ref-dir` (reference taxonomy for specificity; NaN without it),
`--out-dir`, `--results-dir`, `--mock-name`, `--regions`. Rules and truth:
[scoring.md](scoring.md).

## Step 05: real-world tier

`05_realworld_classify.py` runs the classifiers on each study's ASV FASTA and then
checks that every output covers exactly the query set. `05_realworld_evaluate.py`
aggregates ASV counts and the paired shotgun OGU table per sample and rank and
compares them. See [realworld.md](realworld.md).

## Step 06: runtime

`06_runtime.py` samples 100, 500, 1,000, 2,000 and 5,000 sequences, runs each
classifier `--replicates` times and records wall time and peak memory (`psutil`).
Output: `results/runtime/`.

## Parallel runs

Tools receive their own thread count (`--threads`). `03_assign.py --n-jobs K` runs K
commands at once. Separate mock and region cells are independent and can be submitted
as separate batch jobs. QIIME 2 query artifacts are written per classifier
(`rep_seqs_<mock>.<classifier>.qza`) so concurrent jobs do not overwrite each other,
and an artifact older than its FASTA is imported again.

## Re-running part of the matrix

```bash
rm -rf results/mockrobiota-16S/V4/blast
python scripts/03_assign.py --mock-name mockrobiota-16S --regions V4 --classifiers blast
python scripts/04_evaluate.py --mock-name mockrobiota-16S --regions V4
```

Step 02 takes `--force` to replace generated mock data; step 03 overwrites existing
output directories.
