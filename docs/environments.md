# Environments

Scoring and tests need only `requirements-eval.txt`. Running the classifiers needs the
conda environments below. Set their prefixes once, in `conda_envs` of
`configs/pipeline_config.json`; the classifier configs refer to them by key.

| Key in `conda_envs` | Environment | Used by |
|---|---|---|
| `qiime2` | `taxobench-qiime2` | QIIME 2 classifiers, IDTAXA, Kraken 2, MAPseq, output parsers, steps 01 to 06 |
| `spingo` | `taxobench-spingo` | SPINGO |
| `taxonova` | `taxobench-qiime2` with Taxonova 26.10.1 or later | default panel |
| `taxonova_paper` | `taxobench-taxonova-paper` | `classifiers_paper.json` and `classifiers_sweep.json` |

If `conda` is not on `PATH` (for example on some HPC compute nodes), set `conda_bin`
in `pipeline_config.json` or `TAXOBENCH_CONDA_BIN`.

## Scoring only

```bash
python3.11 -m venv .venv && . .venv/bin/activate
pip install -r requirements-eval.txt
```

## taxobench-qiime2

```bash
mamba env create -n taxobench-qiime2 \
    --channel-priority flexible \
    --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-linux-conda.yml

mamba install -n taxobench-qiime2 -y -c bioconda -c conda-forge \
    kraken2 mapseq seqkit art biopython joblib pyyaml psutil \
    bioconductor-decipher r-data.table

# QIIME 2 plugins that can be dropped by later installs
mamba install -n taxobench-qiime2 -y \
    -c qiime2/label/r2024.10 -c conda-forge -c bioconda \
    q2-feature-classifier rescript

# Current Taxonova release (default panel)
mamba run -n taxobench-qiime2 pip install "taxonova>=26.10.1"
mamba run -n taxobench-qiime2 taxonova models download all
```

Flexible channel priority is needed; strict priority does not solve.

**BLAST and `libeutils.so`.** The `entrez-direct` package, pulled in by the install
chain, adds a `lib/libeutils.so` that shadows the library BLAST+ needs, and `blastn`
fails with `symbol lookup error: ... undefined symbol: _ZTIN4ncbi16CESearch_RequestE`.
Move the shadowing file aside:

```bash
mv "$CONDA_PREFIX/lib/libeutils.so" "$CONDA_PREFIX/lib/libeutils.so.entrez-direct-bak"
```

## taxobench-spingo

SPINGO's C library requirements conflict with QIIME 2, so it has its own environment.

```bash
mamba create -n taxobench-spingo -y -c bioconda -c conda-forge spingo
```

## taxobench-taxonova-paper (published configuration)

Taxonova 1.0.0, the version used for the paper, with the v1.0.0 checkpoints placed as
listed in `configs/classifiers_paper.json` (each entry gives the download URL, md5 and
sha256; `03_assign.py` checks the md5 before use).

Taxonova 1.0.0 is not on PyPI, and the `v1.0.0` tag cannot be installed with
`pip install git+...@v1.0.0`: its `setup.py` names a launcher, `bin/taxonova`, that the
tag does not contain. `scripts/install_taxonova_paper.sh` checks out the tag, adds the
launcher (a call to `taxonova.main:main`) and installs the package:

```bash
mamba create -n taxobench-taxonova-paper -y -c conda-forge python=3.10 pandas git
mamba run -n taxobench-taxonova-paper bash scripts/install_taxonova_paper.sh
mamba run -n taxobench-taxonova-paper taxonova classification --help
```

The output parser `src/process_taxonova_output.py` runs in the same environment and
needs only pandas.

## Check

```bash
conda activate taxobench-qiime2
for t in qiime kraken2 mapseq seqkit art_illumina blastn taxonova Rscript; do which $t; done
Rscript -e 'cat("DECIPHER:", "DECIPHER" %in% rownames(installed.packages()), "\n")'
python -c "from qiime2.plugins.feature_classifier.methods import makeblastdb; print('q2-feature-classifier ok')"
python -c "from qiime2.plugins.rescript.pipelines import evaluate_fit_classifier; print('rescript ok')"
taxonova models
conda deactivate

conda activate taxobench-spingo && spingo --help 2>&1 | head -3
```
