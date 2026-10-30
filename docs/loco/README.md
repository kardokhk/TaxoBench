# LOCO: leave-one-clade-out evaluation of novel taxa

## Purpose

Classifiers make two distinct errors: misclassification (MC), naming the wrong taxon
when the true one is in the reference, and over-classification (OC), naming a taxon
when the true one is absent from the reference. The protocol follows the
leave-group-out designs of IDTAXA (Murali et al. 2018) and q2-feature-classifier
(Bokulich et al. 2018). The standard mock and real-world tiers cannot measure OC,
because every true taxon in the queries is also in the reference.

## Definitions (Murali et al. 2018)

For a held-out query *q* with true taxonomy *T_q*, classified at rank L with
confidence *c* against threshold *τ*:

| Condition | Reference contains *T_q* at rank L? | Prediction at rank L | Outcome |
|---|---|---|---|
| 1 | yes | wrong taxon, c ≥ τ | MC error |
| 2 | yes | correct taxon, c ≥ τ | correct |
| 3 | yes | no call (c < τ) | under-classification |
| 4 | no  | any taxon, c ≥ τ | OC error |
| 5 | no  | no call (c < τ) | correct abstention |

## Hold-out design (one per region)

For each region (V4, V3-V4, V5-V7, V1-V9, RRN), one stratified hold-out creates novel
queries at every rank at once:

| Rank | Taxa held out (upper bound) | Query sequences |
|---|---|---|
| Species | 200 | all sequences of those species |
| Genus   | 100 | all sequences of those genera |
| Family  | 50 | all sequences of those families |
| Order   | 25 | all sequences of those orders |
| Class   | 10 | all sequences of those classes |

Counts may be lower for sparse regions. Sampling uses a fixed seed; only taxa with at
least four sequences in the parent rank are eligible, so sister taxa remain.

A query whose family was held out is novel at family, genus and species. At order and
class, where its true taxon is still in the retained reference, it contributes to the
MC denominator. The OC denominator at family is the set of queries whose true family
was held out, and likewise at the other ranks.

## Pipeline

```
scripts/loco/01_make_holdout.py         held-out reference, queries and truth
scripts/loco/02_train_classifier.py     retrain or reindex one classifier on the retained reference
scripts/loco/03_classify.py             classify the queries with one classifier
scripts/loco/04_evaluate.py             MC and OC rates over a threshold sweep
scripts/loco/04c_evaluate_opengates.py  Taxonova curve from open-gate output, one threshold at every rank
scripts/loco/05_plot.py                 curves
```

The default root is `data/loco/` and `results/loco/`. Replicates use their own root
and tag so they never overwrite the default hold-out: `--loco-root loco2k --tag
seed2026` on every step, and `--fixed-higher-ranks <holdout_manifest.json>` on step 01
to redraw only the species stratum. `01_make_holdout.py` refuses to overwrite a
non-empty hold-out.

The eight classifiers of the main panel are supported (BLAST and VSEARCH share one
index, `blast-vsearch`, in step 02). `03_classify.py --maxrejects` sets the VSEARCH
`--p-maxrejects` value; the default is exhaustive search.

## Compute

Per region, approximately: Taxonova training on a GPU 20 to 60
minutes; IDTAXA LearnTaxa about 30 minutes; naive Bayes fit-classifier 1 to 4 hours;
index rebuilds for BLAST, VSEARCH, SPINGO, MAPseq and Kraken 2 5 to 60 minutes. Each
(region, classifier) pair is an independent job.

## References

- Murali A, Bhargava A, Wright ES. *IDTAXA*. Microbiome 2018;6:140. PMID 30092815.
- Bokulich NA et al. Optimizing taxonomic classification of marker-gene amplicon
  sequences with QIIME 2's q2-feature-classifier plugin. Microbiome 2018;6:90.
  PMID 29773078.
- Edgar RC. SINTAX, a simple non-Bayesian taxonomy classifier. bioRxiv 2016.
