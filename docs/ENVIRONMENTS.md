# Conda environments

Two envs total. Recreate from scratch like this:

## 1. `taxobench-qiime2`  (covers 7 of 8 classifiers + tooling)

```bash
# Miniforge (skip if you already have one)
curl -L -o miniforge.sh \
    https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash miniforge.sh -b -p $HOME/miniforge3
source $HOME/miniforge3/etc/profile.d/conda.sh

# QIIME2 amplicon distribution (must use --channel-priority flexible — strict
# priority is unsolvable due to deblur/sortmerna)
mamba env create -n taxobench-qiime2 \
    --channel-priority flexible \
    --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-linux-conda.yml

# Classifier extras
mamba install -n taxobench-qiime2 -y -c bioconda -c conda-forge \
    kraken2 mapseq seqkit art biopython joblib pyyaml psutil \
    bioconductor-decipher r-data.table

# Python deps for taxonova
conda activate taxobench-qiime2
pip install torch tqdm setproctitle

# Install taxonova package (point at the latest taxonova source tree).
# tx_wip_01 is the current default — it ships the comment-2.9 fixes:
#   (a) species-group multi-species labels for region-tied species
#   (b) bf16-AMP / multi-worker / persistent-worker training speedups
#       (toggleable via TAXONOVA_AMP, TAXONOVA_NUM_WORKERS env vars)
#   (c) hierarchical confidence monotonicity at inference (parent ≥ child)
# Install editable so future patches propagate without re-pip-install:
pip install -e /path/to/workdir/taxonova/tx_wip_01/

# QIIME2 plugins that are dropped during transitive reinstalls — restore explicitly
mamba install -n taxobench-qiime2 -y \
    -c qiime2/label/r2024.10 -c conda-forge -c bioconda \
    q2-feature-classifier rescript
```

### Known gotcha — `entrez-direct` clobbers `libeutils.so`

The `entrez-direct` package (pulled in transitively by something in the install
chain) installs `lib/libeutils-1.0.so` and a `libeutils.so` symlink that
**shadows** the `lib/ncbi-blast+/libeutils.so` that BLAST+ needs. `blastn` then
fails with:

```
blastn: symbol lookup error: .../libxobjedit.so: undefined symbol: _ZTIN4ncbi16CESearch_RequestE
```

Fix (already applied in this env — re-apply if you rebuild):

```bash
mv $HOME/miniforge3/envs/taxobench-qiime2/lib/libeutils.so \
   $HOME/miniforge3/envs/taxobench-qiime2/lib/libeutils.so.entrez-direct-bak
```

## 2. `taxobench-spingo`  (only spingo)

SPINGO has C library expectations that conflict with the QIIME2 env, so it
lives on its own.

```bash
mamba create -n taxobench-spingo -y -c bioconda -c conda-forge spingo
```

## Verify

```bash
source $HOME/miniforge3/etc/profile.d/conda.sh
conda activate taxobench-qiime2
for t in qiime kraken2 mapseq seqkit art_illumina blastn taxonova Rscript; do which $t; done
Rscript -e 'cat("DECIPHER:", "DECIPHER" %in% rownames(installed.packages()), "\n")'
python -c "from qiime2.plugins.feature_classifier.methods import makeblastdb; print('q2-fc ok')"
python -c "from qiime2.plugins.rescript.pipelines import evaluate_fit_classifier; print('rescript ok')"
conda deactivate

conda activate taxobench-spingo
which spingo && spingo --help 2>&1 | head -3
```

If all of those print sensible output, you're good.
