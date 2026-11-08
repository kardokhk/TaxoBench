# Scoring rules and ground truth

## Rules

`scripts/04_evaluate.py --rules` selects how calls are scored. The same rules apply
to every classifier; no rule looks at which tool produced a label. The code is in
`src/scoring_v2.py` (v2) and `scripts/04_evaluate.py` (v1).

| `--rules` | abstention (`...; g__X; s__`, `Unassigned`, `unclassified_*`) | group label `s__[A\|B]` | taxon identity in taxon-set metrics |
|---|---|---|---|
| `v2` (default) | no call: a false negative in P/R/F; left out of detection, sensitivity, specificity, Jaccard and L1/L2; reported as `unassigned_fraction` | correct if the truth is a member, at every rank where both are named | rank prefix and name |
| `v2-strict` | as `v2` | a multi-member group is no call at that rank and below | as `v2` |
| `v1` | counted as a taxon of its own, so it is a false-positive taxon in specificity, detection and L1/L2 | credited rank by rank (same outcome as `v2` when the truth has one species) | whole label string |
| `v1-patched` | as `v1`, but dropped from the taxon-set and L1/L2 metrics | as `v1` | as `v1` |

`v1` is the evaluator used for the Taxonova paper. It has two known defects that `v2`
corrects:

- An abstention lowers specificity, so a tool that abstains more can appear to call
  more absent taxa.
- A correct call whose upper ranks are blank (SPINGO writes `d__; p__; ...; g__X;
  s__X y`) never equals the truth string, so it counts as a false-positive taxon and
  a missed taxon.

Under `v2` a tool that abstains often can reach specificity 1.0 with low recall.
Report specificity together with recall or F-measure, not as a ranking on its own.

Under `v2` and `v2-strict` the species rows also carry `n_named_calls`,
`n_group_calls`, `group_share`, `group_size_mean`, `group_size_max`,
`group_truth_in`, `group_truth_out`, `strict_species_f`, `jaccard`,
`unassigned_fraction` and `tp`, `fp`, `fn`.

`scripts/04b_evaluate_sweep.py` (confidence sweeps) and `scripts/loco/04_evaluate.py`
use the v1 rules.

## Ground truth and N bases

The truth builder (`src/process_mock.py`, step 02) treats an `N` in a source sequence
as a mismatch. Up to the published runs it was a wildcard, so a source sequence with
a long run of `N` matched ASVs it does not represent. Pass `--n-wildcard` to
`02_generate.py` or `02b_generate_v5v7.py` (or set `process_mock.N_WILDCARD = True`)
to restore the old rule. `src/rebuild_truth.py` rebuilds one truth file without
QIIME 2:

```bash
python src/rebuild_truth.py \
  --observed data/mock_community/gutmock/V5-V7/paired/rep_seqs_gutmock.fasta \
  --expected data/mock_community/input_data/gutmock_seqs.fasta \
  --taxonomy data/mock_community/input_data/gutmock_taxa.txt \
  --out gutmock_V5-V7_truth.tsv            # add --n-wildcard for the v1 rule
```

With N as a mismatch, labels change in five of the published cells: gut mock V5-V7
(125 of 138 ASVs; 124 deepen from genus to species), *Bilophila* genus to
*B. intracellularis* in biome_digestive_system and biome_sputum at V4 and V3-V4
(one ASV each), and one MBARC-26 whole-16S ASV that becomes `other`. Cells whose truth
is written by the generator rather than by sequence matching (RRN, and whole 16S of
the biome mocks, Mockrobiota-16S, gut mock and vaginal mock) are not affected.

The corrected gut mock V5-V7 truth is
`tests/fixtures/gutmock_V5-V7/data/mock_community/gutmock/V5-V7/trueish_taxonomy_gutmock.tsv`.
The paper scored that cell on the corrected truth; the copy generated with the
wildcard rule is the uncorrected one.

## Reproducing the published values

Use the v1 rules with the v1 truth. The v1 truth treats `N` as a wildcard, so build it
with `02_generate.py --n-wildcard` and `02b_generate_v5v7.py --n-wildcard`; then
replace the gut mock V5-V7 truth with the corrected file above. Without
`--n-wildcard`, three further cells differ from the v1 truth by one label each
(biome_sputum V4, biome_digestive_system V3-V4 and MBARC-26 V1-V9 in the comparison
of the published truth with a rebuild under the N-as-mismatch rule), and scores
under `--rules v1` no longer match the published tables (for example 59 of 440
values in biome_sputum V4).

```bash
python scripts/04_evaluate.py --rules v1 \
  --results-dir <results root> --truth-dir <data root holding the v1 truth> \
  --ref-dir data/ref_dbs --out-dir results/evaluation_v1
```

On 2026-10-31 this rescoring matched all 47 published cells (eight classifiers,
15,040 values) with a difference of 0 at every value, by string comparison after a
stable sort. `tests/test_regression_v1.py` checks the 240 published values of gut
mock V5-V7 on a small fixture in this repository. The v1 code sums the L1 and L2
terms over a Python set, so the sixth decimal of a value can depend on the string
hash seed (seen once in nine runs, one L2 value, 0.070313 for 0.070312); the
published tables print four decimals and are not affected. Set `PYTHONHASHSEED=0`
for identical output across runs.

A cell regenerated from
`mocks/inputs/` need not match the published one: the read simulation settings of
each published cell (`--count`, `--length`, `--paired`) are not recorded in this
repository, and the reads and ASVs also depend on the ART, QIIME 2 and DADA2 versions. The long-read biome cells
cannot be regenerated byte for byte either: the generator used for the paper seeded
its random draws from Python's per-process string hash, and
`02d_generate_longread_biomes.py` uses a stable seed instead.

The Taxonova outputs themselves come from `configs/classifiers_paper.json`
(Taxonova 1.0.0, v1.0.0 checkpoints, 0.5 at every rank above species, 0.7 at
species; install with `scripts/install_taxonova_paper.sh`).

`configs/classifiers_default.json` runs the current Taxonova release instead, whose
outputs differ from the published ones.
