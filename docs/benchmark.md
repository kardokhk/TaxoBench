# Benchmark design

Every classifier runs on the same reads, against references built from the same
source databases, and is scored by the same evaluator against the same truth. The
benchmark has three tiers of community realism and a separate test of
over-classification.

## Amplicon regions

Primers and length windows are in `configs/regions.json`.

| Region | Primers | Length (bp) | Reference |
|---|---|---|---|
| V4 | 515F / 806R (`GTGCCAGCMGCCGCGGTAA` / `GGACTACNVGGGTWTCTAAT`) | 200 to 350 | Greengenes2 |
| V3-V4 | 341F / 805R (`CCTACGGGNGGCWGCAG` / `GACTACHVGGGTATCTAATCC`) | 400 to 500 | Greengenes2 |
| V5-V7 | 799F / 1193R (`AACMGGATTAGATACCCKG` / `ACGTCATCCCCACCTTCC`) | 300 to 500 | Greengenes2 |
| V1-V9 | 27F / 1492R (`AGAGTTTGATCMTGGCTCAG` / `TACGGYTACCTTGTTACGACTT`) | 1,300 to 1,550 | Greengenes2 |
| RRN | none (full-length operon sequences) | 4,600 to 5,300 | GROND (GTDB taxonomy) |

## Tiers

**Tier 1, physical mock communities.** Mixtures with a known species list:
Mockrobiota-16S (concatenated Mockrobiota 16S mocks), MBARC-26 (Singer et al. 2016;
23 bacteria and 3 archaea), a gut mock and a vaginal mock, each at V4, V3-V4, V5-V7
and V1-V9.

**Tier 2, biome-simulated mocks.** For nine environments (digestive system, lung,
sputum, skin, soil, plants, marine, ocean and a mixed control), Greengenes2 reference
sequences are sampled at genus abundances taken from MGnify
(`mgnify/results/2025-08-20_snapshot/`, retrieved with the scripts in `mgnify/`).
Each `biome_*` mock has a full-operon counterpart (`gtdb_*`) for RRN, built from GROND.

**Tier 3, paired amplicon and shotgun studies.** Real samples have no curated truth,
so each amplicon profile is compared with a shotgun profile of the same sample
(Woltka OGU tables). Studies are defined in the `real_world` block of
`configs/mocks.json`:

| Study | Environment | Region | Amplicon | Shotgun |
|---|---|---|---|---|
| Zuo et al. | gut, pediatric ulcerative colitis | V4 | PRJNA759642 | PRJNA759642 |
| Mantri et al. | soil | V3-V4 | PRJNA717918 | PRJNA717813 |
| Virtanen et al. | vaginal | V3-V4 | PRJEB25778 | PRJEB43447 |
| Cao et al. | coral reef | V4 | PRJNA910824 | PRJNA910824 |
| EMP500 forest soil, plant phyllosphere | soil, plant | V4 | PRJEB42019 | PRJEB42019 |

Downloading, trimming, denoising and shotgun profiling are done once per study
([realworld.md](realworld.md), [realworld_amplicon_qc.md](realworld_amplicon_qc.md)).
Reclassifying the ASVs with a new model and comparing again is cheap.

## Ground truth

Short-read cells (V4, V3-V4, V5-V7) are simulated with ART and denoised with DADA2;
each ASV is then matched to its source sequences to give its "trueish" taxonomy.
Long-read cells (V1-V9, RRN) are built by mutating full-length sequences directly, so
their truth is known by construction (`02c_generate_v1v9_direct.py`,
`02d_generate_longread_biomes.py`). See [scoring.md](scoring.md) for how N bases are
treated when matching.

## Metrics

Mock tiers (`scripts/04_evaluate.py`, `src/evaluate_results.py`, `src/scoring_v2.py`),
per rank from class to species:

- precision, recall and F-measure over sequences;
- detection rate (sensitivity): the fraction of true taxa recovered;
- specificity: how often a tool refrains from naming taxa that are absent from the
  truth, over the taxa of the region's reference;
- L1 and L2 error between the true and predicted relative abundances.

Real-world tier (`scripts/05_realworld_evaluate.py`), per sample and rank: Jaccard
index, Bray-Curtis dissimilarity, Aitchison distance (Euclidean distance after a
centred log-ratio transform), precision, recall, F1, L1 and L2.

How abstentions and group labels such as `s__[A|B]` are scored depends on the rule
set; see [scoring.md](scoring.md).

## Over-classification (LOCO)

Standard mocks cannot measure over-classification, because every true taxon is also in
the reference. The leave-one-clade-out experiment removes whole clades from the
reference and measures how often a classifier names a taxon that is absent
(over-classification) against how often it names the wrong one of the taxa present
(misclassification), as a function of the confidence threshold. Protocol:
[loco/README.md](loco/README.md); code: `scripts/loco/`, `src/loco/`.

## Runtime

`scripts/06_runtime.py` samples 100 to 5,000 reads from a mock, runs each classifier
several times and records wall time and peak memory.
