# Real-world tier: paired amplicon and shotgun studies

Real samples have no curated truth. For each study, the amplicon profile of a sample is
compared with the shotgun profile of the same biospecimen. Classification and
comparison are done by two scripts:

- `scripts/05_realworld_classify.py` runs the classifiers on the study's ASVs;
- `scripts/05_realworld_evaluate.py` compares each amplicon profile with its shotgun
  profile.

Preprocessing (download, primer removal, denoising, shotgun profiling) is done once per
study. Reclassifying the ASVs with a new model and comparing again takes minutes.

## Files per study

Paths are relative to `real_world_input_dir` (`data/real_world_input/` by default):

```
rep_seqs_<study>.fasta                          ASV sequences (classifier input)
<study>/dada2_output/asv_table_ids.csv          sample x ASV count table
<study>/woltka_output/ogu_table.biom            shotgun OGU x sample count table
<study>/qiime2_output/exported_taxonomy/taxonomy.tsv   OGU to Greengenes2 taxonomy
<study>/<study>_study/paired/paired_samples.csv amplicon run to shotgun run
```

`paired_samples.csv` needs at least the columns `sample_id`, `amplicon_16s_run` and
`shotgun_wgs_run`. These files, the classifier outputs and the preprocessing records of
the published studies are in the Zenodo deposit ([doi:10.5281/zenodo.23808921](https://doi.org/10.5281/zenodo.23808921)); the accessions of the raw
data are listed in [DATA.md](../DATA.md).

## Registering a study

Add an entry to the `real_world` block of `configs/mocks.json`:

```json
"virtanen": {
    "fasta":          "rep_seqs_virtanen.fasta",
    "study":          "Virtanen et al",
    "region":         "V3-V4",
    "primers":        "341F/785R",
    "amplicon_table": "virtanen/dada2_output/asv_table_ids.csv",
    "wgs_ogu_biom":   "virtanen/woltka_output/ogu_table.biom",
    "wgs_taxonomy":   "virtanen/qiime2_output/exported_taxonomy/taxonomy.tsv",
    "paired_samples": "virtanen/virtanen_study/paired/paired_samples.csv"
}
```

If any of the last four fields is `null`, the study is classified but not compared.

## Preprocessing a new study

The per-study shell scripts used for the paper belong to the archived benchmark data
and are not yet in this repository. The steps, with the tools used, are:

1. **Download** the paired 16S and shotgun runs from ENA or SRA (`prefetch` and
   `fasterq-dump`, or the ENA FTP links) and write `paired_samples.csv`.
2. **Amplicon reads to ASVs.** Remove primers (cutadapt or Trim Galore), then filter,
   denoise and merge with DADA2. Choose `truncLen` so that the read pair still overlaps
   by at least 12 nt, the `mergePairs` default; when it cannot, use the forward read
   only. [realworld_amplicon_qc.md](realworld_amplicon_qc.md) gives the per-study
   choices and the overlap arithmetic. Write the ASVs as `rep_seqs_<study>.fasta` and
   the count table as `asv_table_ids.csv`.
3. **Shotgun reads to an OGU table.** Trim adapters, remove host reads with bowtie2
   against the host genome, align to the Web of Life 2 reference with bowtie2, and
   classify with `woltka classify` into `ogu_table.biom`.
4. **OGU taxonomy.** With QIIME 2, keep the OGUs present in Greengenes2 and export
   their taxonomy to `qiime2_output/exported_taxonomy/taxonomy.tsv`.

Tool versions and filter settings are given in the Methods of the Taxonova paper.

## Running

```bash
python scripts/05_realworld_classify.py --studies zuo mantri virtanen \
    --threads 16 --n-jobs 2
python scripts/05_realworld_evaluate.py --studies zuo mantri virtanen
```

Predictions are written to
`results/real_world/<study>/<classifier>/<parameters>/taxonomy_final.tsv`.
`05_realworld_classify.py` checks that every output covers exactly the query ASVs.
To compare a new Taxonova model, rerun both scripts with `--classifiers taxonova`;
the outputs of the other classifiers are kept.

Outputs of `05_realworld_evaluate.py`, under `results/evaluation/realworld/`:

| File | Content |
|---|---|
| `all_metrics.csv` | one row per study, classifier, parameter set, rank and sample pair |
| `by_metric/<metric>.csv` | mean over sample pairs; rows are study, classifier and parameter set, columns are ranks |
| `per_sample/<study>_<classifier>_<parameters>.csv` | species-rank detail per sample |

## Metrics

Amplicon and shotgun counts are aggregated to each rank and normalized to relative
abundance per sample. Per sample pair and rank:

| Metric | Definition |
|---|---|
| `jaccard` | taxon-set overlap |
| `bray_curtis` | Bray-Curtis dissimilarity of relative abundances |
| `aitchison` | Euclidean distance of centered log-ratio transforms (zeros replaced by 1e-6) |
| `precision`, `recall`, `f-score` | taxa found by amplicon sequencing, with the shotgun taxa as truth |
| `l1_error`, `l2_error` | sum of absolute and of squared abundance differences |
| `weighted_unifrac` | species rank only, with `--unifrac` (needs the Greengenes2 tree; `src/_unifrac_metric.py`) |

A group label such as `s__[A|B|C]` is resolved against the shotgun taxa: its abundance
is split equally among the members found by shotgun sequencing, and is kept as one
taxon when none is found.

## Registered studies

The paper reports Zuo et al. (gut, V4), Mantri et al. (soil, V3-V4) and Virtanen et al.
(vaginal, V3-V4). Three further studies were used during development: Cao et al.
(coral reef, V4) and two EMP500 sub-studies (forest soil and plant phyllosphere, V4).
Accessions are in [DATA.md](../DATA.md).
