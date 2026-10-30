# Real-world paired 16S + shotgun-WGS evaluation

This doc describes (a) the recipe used to onboard a new paired study so it
can be evaluated against shotgun-WGS taxa, and (b) the framework that
ties it all together: `scripts/05_realworld_classify.py` (classification)
+ `scripts/05_realworld_evaluate.py` (WGS comparison).

The whole system is designed so that the *expensive* steps (downloading
raw reads, trimming, DADA2, bowtie2/woltka) are run **once** per study,
their compact outputs are kept in the benchmark data tree (they are not part
of this git repository; see DATA.md), and the cheap steps
(re-classifying ASVs with a refreshed Taxonova model and re-comparing
to WGS) can be re-run any time.

---

## 1. What lives where, per study

```
data/real_world_input/
├── rep_seqs_<study>.fasta        # ASV sequences (DADA2 output)  (kept)
├── rep_seqs_<study>.qza          # same, QIIME2-imported  (kept)
└── <study>/                      # per-study dir
    ├── data/                     # raw reads + reference indexes  (deleted after processing)
    ├── trimmed_16S/              # trim_galore output  (deleted after processing)
    ├── alignments/               # bowtie2 .sam/.bam  (deleted after processing)
    ├── dada2_output/
    │   ├── asv_sequences.fasta
    │   ├── asv_table_ids.csv     # sample × ASV count matrix  (kept)
    │   ├── asv_table.csv         # sample × ASV-sequence count matrix  (kept)
    │   └── filtered/             # intermediate filtered FASTQs  (deleted)
    ├── woltka_output/
    │   └── ogu_table.biom        # OGU × sample count matrix  (kept)
    ├── qiime2_output/
    │   ├── ogu_table.qza
    │   ├── ogu_gg2_filtered.qza
    │   ├── taxonomy.qza
    │   └── exported_taxonomy/taxonomy.tsv  # OGU → GG2 taxonomy lookup  (kept)
    ├── metadata/
    │   └── sample_metadata.tsv
    └── <study>_study/
        ├── data_info/metadata/<bioproject>_runinfo.csv
        ├── paired/paired_samples.csv         # `sample_id, amplicon_16s_run, shotgun_wgs_run, …`
        └── scripts/                          # per-study processing scripts (data tree, not this repository)
            ├── download_sra/
            ├── 16S/                          # trim_galore + DADA2
            └── shotgun_woltka/                # bowtie2 + woltka + qiime2 (GG2)
```

The `paired_samples.csv` file is the key, it's how
`05_realworld_evaluate.py` joins amplicon → shotgun samples
(typically the same biospecimens were sequenced twice).

`configs/mocks.json` `real_world.<study>` registers the paths
(relative to `data_dir/real_world_input/`):

```json
"virtanen": {
    "fasta":           "rep_seqs_virtanen.fasta",
    "study":           "Virtanen et al",
    "region":          "V3-V4",
    "primers":         "341F/785R",
    "amplicon_table":  "virtanen/dada2_output/asv_table_ids.csv",
    "wgs_ogu_biom":    "virtanen/woltka_output/ogu_table.biom",
    "wgs_taxonomy":    "virtanen/qiime2_output/exported_taxonomy/taxonomy.tsv",
    "paired_samples":  "virtanen/virtanen_study/paired/paired_samples.csv"
}
```

If any of `amplicon_table` / `wgs_ogu_biom` / `wgs_taxonomy` /
`paired_samples` is null, evaluation is skipped for that study and only
classification runs.

---

## 2. Onboarding a new paired study

Below is the canonical recipe used for `virtanen`. Replace the BioProject
IDs and run-name patterns to onboard a new study.

### 2.1. Download paired runs from ENA/SRA

```bash
# Pick a working dir on a node with disk headroom (raw data is heavy).
WORKDIR=data/real_world_input/<study>   # relative to the TaxoBench checkout
mkdir -p $WORKDIR/data/{16S_ITS_combined,shotgun} $WORKDIR/metadata
cp <bioproject>_runinfo.csv $WORKDIR/data/metadata/

# Use prefetch + fasterq-dump (or wget on the ENA FTP URLs).
# Save 16S runs to data/16S_ITS_combined/, shotgun to data/shotgun/.
```

Document which 16S run pairs to which shotgun run in
`<study>_study/paired/paired_samples.csv` with at least the columns
`sample_id, amplicon_16s_run, shotgun_wgs_run`.

### 2.2. 16S → ASVs (Trim Galore + DADA2)

```bash
mamba activate trimgalore
bash <study>_study/scripts/16S/01_trim_galore.sh   # parallel; ~2-3 h for 20-30 samples

mamba activate dada2_working
bash <study>_study/scripts/16S/02_dada2_analysis.sh
# tune truncLen for region: V4 ~ c(150, 140); V3-V4 ~ c(270, 200); V5-V7 ~ c(220, 200)
```

DADA2 should drop these into `dada2_output/`:

- `asv_sequences.fasta` (one record per ASV, headers `ASV_1, ASV_2, …`)
- `asv_table_ids.csv`   (sample × ASV-id count matrix)
- `asv_table.csv`       (sample × ASV-seq count matrix; redundant)

### 2.3. Shotgun → OGU table (bowtie2 + woltka against WoL/GG2)

```bash
mamba activate woltka
bash <study>_study/scripts/shotgun_woltka/04_align_reads_woltka.sh   # bowtie2 → BAM
bash <study>_study/scripts/shotgun_woltka/05_run_woltka_classify.sh  # woltka classify → ogu_table.biom
bash <study>_study/scripts/shotgun_woltka/06_qiime2_greengenes2_process.sh  # filter to GG2 + export taxonomy
```

Expected outputs:

- `woltka_output/ogu_table.biom`            (OGU × sample count matrix)
- `qiime2_output/ogu_gg2_filtered.qza`      (GG2-only OGUs)
- `qiime2_output/taxonomy.qza`              (OGU → taxonomy)
- `qiime2_output/exported_taxonomy/taxonomy.tsv` (TSV form, the file we read)

### 2.4. Promote the rep_seqs file into the framework

```bash
cp <study>/dada2_output/asv_sequences.fasta data/real_world_input/rep_seqs_<study>.fasta
qiime tools import --input-path data/real_world_input/rep_seqs_<study>.fasta \
    --type FeatureData[Sequence] \
    --output-path data/real_world_input/rep_seqs_<study>.qza
```

### 2.5. **Storage cleanup: required**

The repo lives on a 500 GB-budget partition. After step 2.3, **delete**:

```bash
rm -rf $WORKDIR/data/shotgun           # raw shotgun FASTQs (>40 GB typical)
rm -rf $WORKDIR/data/16S_ITS_combined  # raw 16S FASTQs (~1 GB)
rm -rf $WORKDIR/alignments             # bowtie2 BAMs (~5 GB)
rm -rf $WORKDIR/trimmed_16S            # trim_galore output
rm -rf $WORKDIR/dada2_output/filtered  # DADA2 filtered intermediates
```

Keep:

- `dada2_output/asv_sequences.fasta`, `asv_table_ids.csv`, `asv_table.csv`,
  `read_tracking.csv`, `analysis_summary.csv`
- `woltka_output/ogu_table.biom`
- `qiime2_output/ogu_table.qza`, `ogu_gg2_filtered.qza`, `taxonomy.qza`,
  `exported_taxonomy/taxonomy.tsv`
- `metadata/`, `<study>_study/` (scripts + paired_samples.csv), `paper/`

### 2.6. Register in mocks.json

Add a `real_world.<study>` entry as in §1. Include the four pairing
fields if WGS is available; otherwise leave them as `null` and the
evaluator will skip the study.

---

## 3. Running classification + evaluation

### 3.1. Classify

```bash
conda activate taxobench-qiime2
cd TaxoBench
python scripts/05_realworld_classify.py --studies virtanen mantri zuo \
    --classifiers all --threads 16 --n-jobs 2
```

Predictions land at
`results/real_world/<study>/<classifier>/<param>/taxonomy_final.tsv`.

### 3.2. Evaluate against shotgun-WGS

```bash
python scripts/05_realworld_evaluate.py --studies virtanen
```

Outputs:

- `results/evaluation/realworld/all_metrics.csv`: tidy long, **per
  paired sample**: study, classifier, param, level, amplicon_sample,
  wgs_sample, jaccard, bray_curtis, precision, recall, f-score,
  l1_error, l2_error, n_amplicon_taxa, n_wgs_taxa.
- `results/evaluation/realworld/by_metric/<metric>.csv`: wide
  pivot: rows = `(study, classifier, param)`, cols = level, values =
  mean across paired samples.
- `results/evaluation/realworld/per_sample/<study>_<classifier>_<param>.csv`
 , species-level per-sample detail (useful for plots / outlier
  inspection).

### 3.3. Re-running after a Taxonova update

When you refresh a Taxonova model under
`data/ref_dbs/<region>/taxonova/`:

```bash
python scripts/05_realworld_classify.py --studies virtanen --classifiers taxonova
python scripts/05_realworld_evaluate.py
```

Only Taxonova's `taxonomy_final.tsv` is overwritten; the other
classifiers' outputs are untouched and the evaluator picks up whatever
predictions exist. This is the intended "fast loop", no re-DADA2,
no re-shotgun.

---

## 4. Metric semantics

| Metric        | Definition                                         |
|---------------|----------------------------------------------------|
| `jaccard`     | per-sample taxon-set overlap at level L            |
| `bray_curtis` | per-sample weighted abundance dissimilarity        |
| `precision`   | amplicon TP / (TP + FP), with WGS as truth         |
| `recall`      | amplicon TP / (TP + FN), with WGS as truth         |
| `f-score`     | 2 P R / (P + R)                                    |
| `l1_error`    | Σ \|p_amp − p_wgs\| over taxa at level L           |
| `l2_error`    | Σ (p_amp − p_wgs)²                                 |

WGS taxa are derived from `ogu_table.biom × taxonomy.tsv` aggregated to
level L; amplicon taxa from `asv_table_ids.csv × <classifier prediction>`
similarly aggregated. Frequencies are per-sample-normalised before
metric calculation. Species-group labels emitted by Taxonova
(`s__[a|b|c]`) are flattened to their first member to avoid
double-counting at species rank.

---

## 5. Registered studies

All six studies in the `real_world` block of `configs/mocks.json` have their
amplicon table, shotgun OGU table and pairing file registered: `zuo` (V4),
`mantri` (V3-V4), `virtanen` (V3-V4), `coral` (V4), `emp500_soil` (V4) and
`emp500_plant` (V4). Accessions are listed in DATA.md.
