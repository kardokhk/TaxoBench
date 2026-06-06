# LOCO — Leave-One-Clade-Out novel-taxon evaluation

## Why this exists

Reviewer-2 comment 2.2: classifiers should be evaluated on **two distinct error
modes** — *misclassification* (MC) and *over-classification* (OC) — using a
leave-one(-many)-group-out protocol following IDTAXA (Murali 2018, PMC6085705)
and QIIME2 (Bokulich 2018, PMID 29773078). The standard taxobench mock and
real-world experiments cannot test OC because every true taxon in the queries
is also present in the reference DB.

## Definitions (Murali 2018)

For a held-out query *q* with true taxonomy *T_q* classified at rank L with
confidence *c* against threshold *τ*:

| Condition | Reference contains *T_q* at rank L? | Predicted at rank L? | Outcome |
|---|---|---|---|
| 1 | yes | wrong taxon, c ≥ τ | **MC** error (mis-classification) |
| 2 | yes | correct taxon, c ≥ τ | correct |
| 3 | yes | no call (c < τ) | under-classification (acceptable) |
| 4 | no  | any taxon, c ≥ τ | **OC** error (over-classification) |
| 5 | no  | no call (c < τ) | correctly abstained |

## Holdout design (one per region)

For each region (V4, V3-V4, V5-V7, V1-V9, RRN), build a single stratified
holdout that creates novel queries at **every rank simultaneously**:

| Rank | Target n taxa held out | Query sequences |
|---|---|---|
| Species | 200 | all sequences of those species |
| Genus   | 100 | all sequences of those genera |
| Family  | 50 | all sequences of those families |
| Order   | 25 | all sequences of those orders |
| Class   | 10 | all sequences of those classes |

Caps are upper bounds — actual counts may be lower for sparse regions.
Sampling is reproducible: seed fixed; only taxa with **≥ 4 sequences in the
parent rank** are eligible (so removing them still leaves enough sister taxa).

A held-out sequence *q* whose family was held out is novel at family,
genus, and species levels (its order/class are still present). At order/class
levels (where its true taxon IS present in retained DB) *q* contributes to the
**MC denominator** at those ranks. The OC denominator at family is sequences
whose true family was held out, etc.

## Pipeline

```
01_make_holdout.py     # generate held-out reference + query/ground-truth files
02_train_classifiers.py # retrain/reindex all 8 classifiers on held-out ref
03_classify.py         # classify queries with all 8 classifiers
04_evaluate.py         # MC + OC rates at threshold sweep + curves
```

Outputs land under `results/loco/<region>/`.

## Compute notes

Per region: ~20–60 min Taxonova GPU training; idtaxa LearnTaxa ~30 min CPU;
naive-Bayes fit-classifier ~1–4 h CPU; BLAST/VSEARCH/SPINGO/MAPseq/Kraken2
index rebuilds 5–60 min CPU. SLURM array jobs parallelize across regions.

## References

- Murali A, Bhargava A, Wright ES. *IDTAXA*. Microbiome 2018;6:140. PMID 30092815.
- Bokulich NA et al. *Optimizing taxonomic classification of marker-gene
  amplicon sequences with QIIME 2's q2-feature-classifier plugin*. Microbiome
  2018;6:90. PMID 29773078.
- Edgar RC. *SINTAX, a simple non-Bayesian taxonomy classifier*. bioRxiv 2016.
