# Development log

Kept while the benchmark was built for the Taxonova paper (April to October 2026).
It is a historical record: paths, script names and review-round references in it
refer to the development tree and may no longer exist. Current changes are in
[CHANGELOG.md](../CHANGELOG.md).

## 2026-10-29: scoring v2 (default), N-as-mismatch truth, path and SPINGO fixes

* **Scoring rules** (`scripts/04_evaluate.py --rules`, `src/scoring_v2.py`).
  New default `v2`: an abstention is no call (false negative in P/R/F, never a
  taxon in detection, sensitivity, specificity, Jaccard or L1/L2; share reported
  as `unassigned_fraction`); a group label is correct iff the truth is a member
  at every rank where both are named; taxon identity in the taxon-set metrics
  is rank prefix + name. `v2-strict` treats a multi-member group as no call.
  `v1` keeps the published code path unchanged and reproduces the published
  values; `v1-patched` is v1 with abstentions dropped from the taxon-set
  metrics only. New options `--truth-dir`, `--ref-dir`, `--out-dir`.
* **SPINGO taxon identity.** Under v1 a correct SPINGO call (blank upper ranks)
  never equalled the truth string and counted as a false-positive taxon plus a
  missed taxon. Fixed in v2 by comparing rank + name.
* **Specificity sentinel.** A missing reference universe now gives NaN rather
  than an empty string, which made the `by_metric` pivot fail.
* **Truth builder.** `src/process_mock.py` treats N in a source sequence as a
  mismatch (was a wildcard). Matching code moved to `src/truth_match.py`
  (no QIIME 2 dependency); `src/rebuild_truth.py` rebuilds one truth file.
  `--n-wildcard` on `02_generate.py` / `02b_generate_v5v7.py` restores the old
  rule. Labels change in 5 published cells, among them gutmock V5-V7 (125 of
  138 ASVs, the correction already used for the paper's round-3 tables).
* **Paths.** `configs/pipeline_config.json` uses relative paths, resolved
  against the repo root; `TAXOBENCH_<KEY>` environment variables override a key.
* **Tests** (`tests/`, 24 tests, one of which needs the data tree): abstention
  formats of every tool, group in/out, strict, SPINGO format, N rule, fast
  search equal to brute force, gutmock V5-V7 truth rebuild, and v1 regression
  against the 240 published gutmock V5-V7 values.

## 2026-05-10 — add MBARC-26 to v_mocks panel; reject CAMI2; drop ZymoBIOMICS

* **Adds**: `mbarc26` (Singer 2016 Sci. Data, 23 bacteria + 3 archaea,
  16S fetched from NCBI RefSeq NR_ records) to the `v_mocks` matrix.
  Evaluated across V4 / V3-V4 / V5-V7 / V1-V9 (no RRN — no published
  full-operon truth).
* **Rejects** CAMI2: shotgun-only challenge, no curated 16S reference +
  amplicon truth. (MBARC-26 is itself an ingredient of CAMI2, so adding
  CAMI2 alongside would double-count.)
* **Drops** ZymoBIOMICS D6300 (built and evaluated, then removed): the
  chkMocks D6300 16S sub-region sequences match gg2 references labelled
  with sister species at higher identity than the canonical-name gg2
  references (S. aureus → argenteus, B. subtilis → inaquosorum,
  S. enterica → bongori, P. aeruginosa → paraeruginosa, L. monocytogenes →
  marthii, E. coli → Shigella sonnei). All 8 gg2-trained classifiers in
  the panel systematically miscall species — a reference-DB nomenclature
  mismatch on the chkMocks data, not a pipeline issue. Removed all zymo
  data + results + config entries.
* **New code**:
  - `scripts/01b_prepare_extra_mocks.py` — fetches + formats MBARC-26
    `data/mock_community/input_data/mbarc26_{seqs.fasta,taxa.txt}` via
    NCBI eutils.
  - `scripts/02c_generate_v1v9_direct.py` — V1-V9 generator for curated
    mocks (no MGnify, no backbone — pool is the curated species fasta).
    Same structure as the env-biome regenerator but with mutation rate
    set to **0** to match mockrobiota-16S's recipe (verified empirically:
    mockrobiota's V1-V9 env-seqs are 100% identical to input fasta after
    primer extraction). Deterministic seed via `zlib.crc32` (not Python's
    process-randomized `hash()`). Lives in taxobench (not taxonova) so
    the production tx_wip_05 install is untouched per ops policy.
* **Config**: new `extra_mocks` block in `configs/mocks.json`;
  `mbarc26` appended to `matrix.v_mocks`.
* **Taxonomy**: lineages resolved against
  `data/ref_dbs/Greengenes2_backbone/2024.09_taxonomy.tsv` (gg2-canonical
  with bracket-form where multiple gg2 labels apply, e.g. E. coli ↔
  Shigella). Two MBARC species absent from gg2 r2024.09
  (`Salmonella enterica`, `Thermobacillus composti`) get genus-only
  fallback lineages — known gg2 gaps, not mock defects.
* **Pipeline uniformity**: identical recipe across mockrobiota-16S,
  biome_*, and mbarc26. Truth = single gg2-canonical species per ASV
  (mockrobiota convention); the bracket-aware evaluator
  (`src/evaluate_results.py:_expand_group_label`) credits taxonova's
  bracket predictions when the truth species is in the bracket.
* **Constraints honoured**: no edits to taxonova source trees
  (`tx_wip_05/`, `tx_wip_07/`, `tx_wip_08/`); no edits to model checkpoints
  under `data/ref_dbs/`. All work is in `taxobench_final/{scripts,configs,data,docs}/`.
* **Detail write-up**:
  `corrections/comments_addressed/extra_mocks_mbarc26.md`.

## 2026-05-09 — env biome V1-V9 + RRN data fix (200 bp truncation root cause)

* **Root cause**: comment-1.2's add of marine/ocean/plants biomes ran V1-V9 and
  RRN through `scripts/02_generate.py`, which uses ART (Illumina simulator,
  caps at ~300 bp/read). Resulting rep_seqs were 200 bp uniform; the V1-V9 +
  RRN taxonova models, trained on full-length references, collapsed on these
  inputs. taxonova species F dropped to 0.04-0.13 vs other classifiers'
  0.52-0.87 — pure pipeline asymmetry, not a model deficiency. Existing
  biomes (digestive/lung/skin/soil/sputum/mixed) avoided this because their
  V1-V9 rep_seqs were generated 2026-04-26 by the manuscript-documented
  long-read recipe (mutate gg2 backbone at 0.002-0.02), not via ART.
* **Affected cells**: biome_marine/ocean/plants × V1-V9; gtdb_marine/ocean/plants × RRN.
* **Fix**: new `regenerate_env_biome_repseqs.py` (in
  `/path/to/workdir/taxonova/tx_wip_08/scripts/`) implements the
  manuscript recipe — sample up to 500 sequences from gg2 backbone (V1-V9) or
  GROND (RRN) weighted by MGnify genus abundance, generate 10 env-seqs
  replicates with per-base mutation rate uniform in [0.002, 0.02],
  dereplicate via MD5(seq) hash to rep_seqs (DADA2-style hashed IDs).
  Applied to all 6 cells; broken 200 bp legacy files moved aside to
  `<biome>/<region>/single/rep_seqs_<biome>.fasta.broken_200bp_aside_2026-05-09`.
* **End-to-end verification (`scripts/04_evaluate.py`, FINAL — full 8-classifier
  panel including RRN vsearch; all SLURM jobs 9578-9583 completed, gtdb_marine
  RRN vsearch took 2d 15h)**: taxonova species F before → after on the 6 affected
  cells: biome_marine V1-V9 0.042 → 0.973 (#1 of 8; idtaxa 0.756, vsearch 0.697);
  biome_ocean V1-V9 0.059 → 0.987 (#1; idtaxa 0.905, vsearch 0.673);
  biome_plants V1-V9 0.073 → 0.965 (#1; idtaxa 0.818, vsearch 0.742);
  gtdb_marine RRN 0.127 → 0.982 (#1; vsearch 0.878, blast 0.863);
  gtdb_ocean RRN 0.129 → 0.990 (#1; blast 0.843, vsearch 0.824);
  gtdb_plants RRN 0.134 → 0.962 (#1; vsearch 0.853, blast 0.846).
  Taxonova goes from catastrophically last to first on every
  cell with the full panel; margin vs next-best 0.07-0.28 F. RRN vsearch
  jobs are slow (~5kb queries × 164k 5kb refs); they were dispatched to
  cpuq-long with --threads 20, --n-jobs 2 via slurm/assign_one_subset.sbatch.
* Pre-eval inference-confidence sanity check on biome_marine V1-V9: genus mean
  conf 0.65 → 0.997, species 0.66 → 0.94. Lengths full ~1490 bp V1-V9 /
  ~4900 bp RRN.
* **Regression guards** (so this can't recur):
  - `scripts/02_generate.py` errors out for `--regions V1-V9` or `RRN`,
    directing callers to the regenerator.
  - `slurm/generate_one.sbatch` errors out for `REGION=V1-V9` or `RRN`.
  - `docs/PIPELINE.md` adds a "Long-read regions" section with the recipe.
  - `corrections/comments_addressed/comment_1_2_biomes_and_realworld.md`
    updated to mark the V1-V9 caveat RESOLVED with the fix recipe.

## 2026-05-05 — real-world evaluator: group-label resolution aligned with mock/biome evaluator

`scripts/05_realworld_evaluate.py` now mirrors the comment-2.9 group-label
semantics of `evaluate_classification` in `src/evaluate_results.py` so
that taxonova's species-group predictions (`<rank>__[a|b|c]` from the
species-group models) are scored consistently across the mock/biome and
real-world panels.

* **New helpers**: `_expand_group_label` (cartesian expansion of bracket
  predictions across any rank — generalises to multi-rank brackets) and
  `_resolve_group_keys` (truth-aware re-keying applied per paired sample
  before any metric is computed).
* **Resolution rule** (one ASV = one prediction unit):
  - bracket prediction `<rank>__[a|b|c]` at frequency `f`;
  - if any of {a,b,c} is in the WGS truth set, `f` is routed to the
    matching member(s) (split equally if more than one matches) — full
    TP credit goes to the truth;
  - if none match, the bracket is kept as a single bucket — counts as
    one false positive, not N.
* **Bracketless predictions** (every other classifier + the non-group
  taxonova model) pass through unchanged. Sanity-checked: with the
  resolver enabled, all seven non-taxonova classifiers produce
  byte-identical numbers vs. running without it.
* **Verification**: 8 unit tests covering bracketless, single-match,
  no-match-keeps-bucket, multi-match-split, higher-rank brackets,
  truncation, mass-conservation, and multi-bracket cases.
* **Impact on taxonova species-level scores** (means across paired
  samples; non-taxonova classifiers unchanged):
  - **zuo** (gut V4, 42 samples): F-score 0.098, BC 0.612, precision
    0.656 — taxonova sweeps the leaderboard.
  - **emp500_soil**: precision 0.263 — taxonova top-1 on precision;
    F-score 0.0494 (#2, edges naive-bayes 0.0493).
  - **coral**: precision 0.239 (was 0.210), F-score 0.032 (was 0.027).
  - **emp500_plant**: precision 0.329, F-score 0.0998.
  - **mantri**: precision 0.384 — taxonova top-1.
  - **virtanen**: F-score 0.027, BC 0.877 — taxonova top-1 on both
    (no group labels appeared in those predictions, so unchanged
    from prior).

## 2026-05-04 — EMP500 environmental studies added (comment 1.2 amendment b)

* **Two new paired-WGS real-world studies** added to broaden the panel
  with under-characterised environmental biomes (reviewer 1.2):
  - **`emp500_soil`** — EMP500 sub-study `myrold.59` (forest soil),
    22 paired biospecimens, 16S V4 single-end MiSeq + Illumina
    NovaSeq shotgun.
  - **`emp500_plant`** — EMP500 sub-study `berry.2` (plant
    phyllosphere / berry surface), 20 paired biospecimens, same
    sequencing combo.
* **Pipelines** added: `slurm/emp500_amplicon_pipeline.sbatch` (parametric
  via `STUDY=emp500_soil|emp500_plant`), `slurm/emp500_wgs_pipeline.sbatch`
  (same), `slurm/emp500_woltka_recover.sbatch` (BAM→SAM recovery).
* **EMP500 amplicon-pick fix**: PRJEB42019 mixes 16S, 18S, and ITS
  amplicons under one `library_strategy=AMPLICON` label. Per-biospecimen
  picks were initially the deepest amplicon run, which sometimes hit
  ITS reads. Fixed by filtering on `library_construction_protocol`
  matching the two strings *"Illumina EMP protocol 515fbc[, ]806r
  amplification of 16S rRNA V4"*. Verified by inspecting ASV_1 in both
  studies (now starts with canonical post-515F `TACAGAGG…`, not the
  ITS1 18S/5.8S junction `AAGTCGTAACAAGGT`).
* **woltka SAM/BAM gotcha documented**: woltka's `infer_align_format`
  reads bytes as text and chokes on BAM's gzip magic. The mantri
  template kept `.sam` files alongside `.bam`; my initial EMP500
  draft deleted SAMs after BAM conversion → woltka failed. Patched
  to skip BAM compression entirely (woltka reads SAM directly, then
  the whole `alignments/` dir is dropped). Recovery script handles
  pre-existing BAMs.
* **Real-world panel now has 6 studies × 8 classifiers × 6 ranks =
  6 096 paired-sample rows** in `results/evaluation/realworld/`.
* **`configs/mocks.json`** `real_world.{emp500_soil,emp500_plant}`
  blocks fully populated.
* **`corrections/comments_addressed/comment_1_2_biomes_and_realworld.md`**
  amended with both new studies in the species-level table and
  per-study leaderboard.

## 2026-05-04 — real-world amplicon paired-end audit (comment 1.2 amendment a)

* **Audited the four real-world amplicon DADA2 pipelines** for
  trimming / filtering / merge correctness; documented in
  `docs/REALWORLD_AMPLICON_QC.md`.
* **Bug fixed: zuo single-end → paired-end.** The original
  `zuo_amplicon_pipeline.sbatch` ran DADA2 single-end with
  `truncLen=110` on the (incorrect) reasoning that 2×150 + 39 bp
  primers leaves no overlap. Zuo reads are *already* primer-trimmed
  by the submitter, so R1=R2=150 bp covers the V4 insert end-to-end
  (overlap = 150+150−252 = **48 bp**). New
  `slurm/zuo_amplicon_pipeline_paired.sbatch` runs paired-end with
  `truncLen=c(145,140)` → **3 218 ASVs at 252–253 bp** (was 1 769 at
  110 bp). Re-classified all 8 classifiers and re-ran
  `05_realworld_evaluate.py`.
* **Coral confirmed correct as single-end.** With primers in the read
  and only 2×150 chemistry, post-cutadapt overlap is 9 bp, below
  DADA2's `minOverlap=12`. An attempted paired-end variant dropped
  every read in `filterAndTrim` (truncLen too long for the
  131/130 bp inputs). The original single-end pipeline is the right
  choice; documented in `docs/REALWORLD_AMPLICON_QC.md` §2.3.
* **Mantri & virtanen verified correct** (V3-V4 paired-end; ASVs at
  ~410 / ~445 bp respectively).
* **`corrections/comments_addressed/comment_1_2_biomes_and_realworld.md`**
  amended at top with the zuo updates; species-level table and
  per-study leaderboard refreshed. **`results/evaluation/realworld/`
  refreshed → 4 080 paired-sample rows.**

## 2026-05-04 — comment 1.2 follow-ups (coral + Zuo paired-WGS)

* **Coral Reef paired-WGS demo** added (PRJNA910824 — Florida Coral
  Reef, 15 paired biospecimens via `16S_<id> ↔ MG_<id>` alias matching).
  V4 region (515FY/806RB), iSeq 2x150, single-end DADA2 (paired merge
  impossible: 252 bp V4 + 39 bp primers = 290 bp ≈ paired read budget).
  4,876 ASVs. New SLURM pipelines `slurm/coral_{amplicon,wgs}_pipeline.sbatch`.
* **Zuo paired-WGS reprocess** (PRJNA759642 — pediatric ulcerative
  colitis vs healthy gut). User's prior Zuo data was lost; reprocessed
  from scratch using the same template. 42 paired biospecimens via
  `metadata.csv` from `github.com/wbb121/IBD-data-analysis`. Reads are
  pre-primer-trimmed by submitter (start with `TCCGTAGG`, not 515F's
  `GTGCCAG`); cutadapt config skips primer removal accordingly.
  1,769 single-end V4 ASVs across 42 samples. New SLURM pipelines
  `slurm/zuo_{amplicon,wgs}_pipeline.sbatch`.
* **Both new studies** survived an OOM/HTML-error class of bugs:
  - bowtie2 OOM on 128 GB nodes → switched all WGS pipelines to
    `cpuq-short` 1 TB nodes with `--mm` shared-memory index (no peak
    memory issue with ~50 GB RSS for the index map).
  - ENA FTP sometimes returns an HTML 404 directory listing as the
    `_2.fastq.gz` for runs that are submitted single-end → added a
    `gzip -t` integrity check after each download (corrupt → fall
    back to single-end).
* **`results/evaluation/realworld/`** now spans **4 studies × 8
  classifiers × 6 ranks = 4,032 paired-sample rows**.
* **`configs/mocks.json`** `real_world.{coral,zuo}` blocks fully
  populated (region, primers, amplicon_table, wgs_ogu_biom,
  wgs_taxonomy, paired_samples).

## 2026-05-03 — comment 1.2 (broaden biomes + add a 2nd paired-WGS demo)

* `scripts/05_realworld.py` → renamed `scripts/05_realworld_classify.py`
  (classification only — same behaviour as before).
* **NEW** `scripts/05_realworld_evaluate.py` — compares each
  `(study × classifier × param)` ASV classification against the paired
  shotgun-WGS taxa table, reports Jaccard + Bray-Curtis (reviewer 1.8) +
  precision/recall/F + L1/L2 per paired-sample, and writes mean-aggregated
  wide CSVs per metric under `results/evaluation/realworld/`. Smoke-tested
  on virtanen (1,008 per-sample rows across 8 classifiers × 6 levels × 21
  paired samples).
* `configs/mocks.json` — `real_world.<study>` entries gained four
  optional fields (`amplicon_table`, `wgs_ogu_biom`, `wgs_taxonomy`,
  `paired_samples`). When all four are populated the evaluator runs;
  otherwise it's silently skipped. Currently wired for `virtanen`.
* `configs/mocks.json` — registered `biome_marine`, `biome_plants`,
  `biome_ocean` (V-region) and `gtdb_marine`, `gtdb_plants`, `gtdb_ocean`
  (RRN counterparts). Input fastas/taxa already exist under
  `data/mock_community/input_data/`; mock generation per region
  (`02_generate.py` + `02b_generate_v5v7.py`) and downstream
  classification/evaluation still need to be run.
* **NEW** `docs/REALWORLD.md` — recipe for onboarding a new paired study
  (download → trim_galore → DADA2 → bowtie2/woltka → register in
  `mocks.json`), storage-cleanup policy (delete shotgun + raw +
  intermediates after woltka completes), and the fast
  re-classification/re-evaluation loop intended for Taxonova model
  refreshes.
* `data/real_world_input/virtanen/` cleanup: removed
  `data/shotgun` (52 GB), `data/16S_ITS_combined` (1.1 GB),
  `alignments/` (4.9 GB), `trimmed_16S/` (615 MB),
  `dada2_output/filtered/` (104 MB) — ~58 GB freed; only ASV table /
  OGU biom / GG2-mapped taxonomy / scripts / metadata retained.

## 2026-05-02 — comment-2.9 fixes (Taxonova species-group + monotonic confidence + speedups)

Three coordinated changes that flip on automatically next time `03_assign.py`
or `04_evaluate.py` is run:

1. **Default Taxonova model = species-group model** (per region). The
   pre-trained models under `data/ref_dbs/<region>/taxonova/` were renamed
   from `gg2_<region>_model.pt` to `gg2_<region>_species_groups_model.pt`.
   The previous models are kept alongside as `gg2_<region>_model.pt.baseline_bak`
   for reproducibility. The `03_assign.py` model glob (`*model.pt`) ignores
   the `.baseline_bak` suffix and so always picks the new model.

2. **Default Taxonova package = `tx_wip_01`** in the `taxobench-qiime2` env
   (was `taxonova_113`). `tx_wip_01` ships:
   - Hierarchical confidence monotonicity at inference: each rank's reported
     confidence is `min(raw_conf, parent_conf)` so `domain ≥ phylum ≥ ... ≥
     species` always holds. Fixes the long-standing `species_conf > genus_conf`
     anomaly (5,544 violations in the prior baseline outputs).
   - bf16 mixed-precision training (`TAXONOVA_AMP=1` default), multi-worker
     DataLoaders (`TAXONOVA_NUM_WORKERS=4` default), persistent workers,
     `set_to_none` zero_grad, non-blocking H2D — ~3-4× per-epoch training
     speedup on L40S GPUs.

3. **Group-aware evaluator** (`src/evaluate_results.py`): the per-sequence
   `precision_recall_fscore` → `evaluate_classification` chain now expands
   `s__[a|b|c]` group labels to a set and counts a TP at species rank iff
   the expected species is a member. Backward-compatible: single-species
   labels behave identically. The wide-format CSVs (`by_metric/*.csv`) and
   the long-format `all_metrics.csv` automatically reflect this.

**Headline V4 mock results** (post-update, species-level F1 mean across
the 7 biome mocks): 0.867 (baseline) → 0.928 (species-group + group-aware
eval). V5-V7: 0.700 → 0.936. V3-V4 / V1-V9: unchanged (already-resolved
regions; species-group fix correctly inactive). Rationale + per-mock
tables in `corrections/comments_addressed/comment_2_9_ood_and_ambiguous_species.md`.

**To revert** (e.g. for an apples-to-apples re-run vs. the original Table 1):
```bash
cd data/ref_dbs
for r in V4 V3-V4 V5-V7 V1-V9; do
    KEY=$(echo $r | tr '[:upper:]-' '[:lower:]_'); KEY=${KEY//_/}
    mv $r/taxonova/gg2_${KEY}_species_groups_model.pt $r/taxonova/gg2_${KEY}_species_groups_model.pt.disabled
    mv $r/taxonova/gg2_${KEY}_model.pt.baseline_bak    $r/taxonova/gg2_${KEY}_model.pt
done
# and for the evaluator: git checkout the prior src/evaluate_results.py.
```

## What changed vs. taxobench_11

### Removed permanently

- **3 classifiers**: `rdp`, `mothur`, `emu` (their `process_*_output.py`,
  `build_*_command()` builders, and ref DB symlinks all gone).
- **ITS / UNITE / fungal** logic — `extract_its.py`, `parse_mmseqs_output_ITS.R`,
  ITSx wrappers, all references to `ITS1`, `ITS2`, `ITS12` regions, and the
  `mockrobiota-ITS` flow.
- **Cross-validation** — every CV branch in `evaluate_results.py`, `process_mock.py`.
- **All matplotlib / seaborn plotting** — outputs are CSV-only.
- **Combined multi-level result tables** — every metric is now its own CSV
  (`{level}_{metric}.csv`).
- ~~L1 / L2 abundance-error metrics~~ — re-added 2026-04-29 (see below).
- **Hardcoded paths** — `/home/kardokhk/...` references replaced with config-driven
  paths under `taxobench_final/`.

### Restructured

- **One canonical layout**: `notebooks/ src/ scripts/ configs/ slurm/ docs/ data/
  results/ logs/`.
- **Standalone CLI runners** under `scripts/` (`01_prepare.py` … `06_runtime.py`),
  each driven by `argparse`. The original notebook `.ipynb` files are NOT
  carried over — the same logic now lives in importable modules + CLI scripts.
- **Two conda envs** (`taxobench-qiime2`, `taxobench-spingo`) instead of three
  loosely-named ones.
- **Reference DBs** physically migrated into `data/ref_dbs/` from
  `taxobench_11/benchmark/data/ref_dbs/` (V4, V3-V4, V1-V9, RRN [from
  `gtdb_rrn`], Greengenes2_backbone [from `Greengens2_backbone`]). ITS/UNITE
  dirs left in `taxobench_11` (out of scope).
- **Mock community data** physically migrated into `data/mock_community/`:
  `mockrobiota-16S`, the 6 in-scope biomes (`biome_{sputum,soil,skin,lung,
  digestive_system,mixed}`) and their `gtdb_*` RRN counterparts, plus
  `input_data`, `gutmock`, `vaginalmock`, `biomes_source`. Out-of-scope
  biomes (fecal, oral, marine, ocean, saliva, vagina, plants, human_*),
  `mockrobiota-ITS`, `prepare_mock_files`, and stale `results/` left in
  `taxobench_11`.

### Added

- 4 config files in `configs/`: `pipeline_config.json`, `classifiers_default.json`,
  `mocks.json`, `regions.json`.
- Slurm submission templates: `slurm/assign_one.sbatch`,
  `slurm/assign_full_matrix.sh`, `slurm/realworld.sbatch`, `slurm/runtime.sbatch`.
- Real-world flow now properly handles 3 studies (Zuo one, Zuo two, Mantri)
  instead of 1.
- Runtime / RAM benchmark module (`06_runtime.py`) wired into the main pipeline
  layout instead of living as a side notebook.
- Docs: `README.md`, `docs/PIPELINE.md`, `docs/ENVIRONMENTS.md`, `docs/CHANGELOG.md`.

## Smoke-test results (Mockrobiota-16S × V4)

All 8 classifiers ran end-to-end and produced 95-line `taxonomy_final.tsv`:

| Classifier  | Status | Notes                                                      |
|-------------|--------|------------------------------------------------------------|
| taxonova    | yes      | Required newer model from `taxobench_20/...gg2_v4_model.pt` (the version-1 model in `taxobench_11/` is incompatible with taxonova_113) |
| vsearch     | yes      | ~30 s                                                      |
| blast       | yes      | Required `libeutils.so` shadow-fix (see ENVIRONMENTS.md)   |
| naive-bayes | yes      |                                                            |
| idtaxa      | yes      |                                                            |
| spingo      | yes      | Runs in `taxobench-spingo` env, processing in main env     |
| kraken2     | yes      |                                                            |
| mapseq      | yes      |                                                            |

Evaluate stage: 20 per-(level × metric) CSVs written under
`results/mockrobiota-16S/V4/`. Sample `Genus_f-score.csv`:

```
method,...,confidence_threshold-0.7,...
blast,...,0.795
idtaxa,...,0.926
kraken2,...,0.561
mapseq,...,0.911
naive-bayes,...,0.951
spingo,...,0.933
taxonova,...,0.936
vsearch,...,0.81
```

Real-world: `zuo_one × V3-V4 × vsearch` → success.
Biome: `biome_sputum × V4 × vsearch` → success.
RRN: `gtdb_sputum × RRN × vsearch` → submitted; takes ~10 min on full-length
GROND seqs — same machinery as V-region runs.

## Things flagged for the user

- **`01_prepare.py`** does NOT invoke `mmseqs2` — for fresh Mockrobiota datasets
  it just writes empty taxa files. The original notebook's `mmseqs2 +
  parse_mmseqs_output_16S.R` flow was not ported because all current mocks are
  already prepared and symlinked. If you ever need to ingest brand-new
  Mockrobiota mocks, port that step from `taxobench_11/scripts/prepare_mock_community.ipynb`
  or call the original notebook once.
- **`02_generate.py`** — same story; current mocks are already-generated. The
  primer-extraction call is defensive (`hasattr(extract_reads, 'extract_amplicon')`)
  with a copy-fallback. For brand-new mocks you may need to align to
  `extract_reads`'s actual public API.
- **Old taxobench_*/ directories** under `/path/to/workdir/taxobench/`
  are now safe to delete. All in-scope data has been physically migrated
  into `taxobench_final/data/` (no remaining symlinks to `taxobench_11`).
  The taxonova package source lives at
  `/path/to/workdir/taxonova/taxonova_113/`.

  *(Historical note: the project root was originally `taxonova_bench/`; it
  was renamed to `taxobench/` on 2026-04-29 to streamline naming. References
  to the old path in older changelog entries below have been preserved
  verbatim, but no current scripts or configs depend on it.)*

## Data migration (2026-04-29)

All symlinks in `data/` replaced with the actual directories via `mv` (same
filesystem, instant). Sources drained from:

- `taxobench_11/benchmark/data/ref_dbs/` → `data/ref_dbs/`
- `taxobench_11/benchmark/data/mock-community/` → `data/mock_community/`
- `taxobench_11/new_data/` → `data/real_world_input/`
- `taxobench_11/benchmark/softy/` → `data/softy/`
- `taxobench_20/taxonova_113/` → repo-root `taxonova_113/`

`taxobench_final/data/` now holds ~31 GB of real data (no symlinks into
sibling taxobench_* dirs).

### L1 / L2 abundance-error metrics re-added (2026-04-29)

Ported back from `taxobench_11/scripts/evaluate_mock_community.ipynb`. New
implementation in `scripts/04_evaluate.py::_abundance_norms`: for each level,
count occurrences of each truncated taxon in true vs. predicted lists,
normalize to relative frequencies, then
`L1 = Σ |p_exp − p_obs|`, `L2 = Σ (p_exp − p_obs)²` over the union of taxa
(L2 kept squared, not sqrt — matches the old convention).

Old version computed L1/L2 only at species level via biom tables. New
version covers **all 5 levels** (Class–Species) using the same TSV pipeline
the other metrics already use — no biom dependency.

Output: `results/evaluation/by_metric/{l1_error,l2_error}.csv` plus
`l1_error`/`l2_error` rows in `all_metrics.csv`.

### MGnify retriever (2026-04-29)

Imported the MGnify top-taxa retrieval helper into `mgnify/`. Used to pull
the most-abundant taxa per rank for a biome — this feeds the manual
construction of new `biome_*` mocks under `data/mock_community/`.

Reorganized from the original ad-hoc layout (7 numbered scripts, 18
single-biome timestamp dirs) into:

- `mgnify/list_biomes.py` — refresh the lineage catalogue (`biomes.csv`).
- `mgnify/fetch_top_taxa.py` — one CLI-driven script (was `2_*`/`5_*`/`6_*`
  + duplicates), argparse-based, function-style, with `--biomes`,
  `--lineage`, `--max-samples`, `--top-n`, `--run-tag`, `--verbose`.
- `mgnify/biomes.csv` — renamed from `mgnify_environments.csv`.
- `mgnify/results/2025-08-20_snapshot/` — old per-biome dirs consolidated
  here (the source for current `biome_*` mocks).
- Removed: 3 debug-only scripts, 1 duplicate of script 6, 18 fragmented
  result dirs.
- See `mgnify/README.md` for usage.

Smoke-tested against the live MGnify v1 API (2 marine samples → 7 rank CSVs
in 17 s).

### Training recipes (also 2026-04-29)

Kraken2 and mapseq DB-training pipelines moved into `training/`:
- `taxobench_13_kraken/kraken2_train_db/` → `training/kraken2/`
- `taxobench_15_mapseq/train_mapseq/`     → `training/mapseq/`

Not integrated into `scripts/03_assign.py` — these are parked recipes for
infrequent DB refreshes. Deployed kraken2/mapseq dbs already live under
`data/ref_dbs/<region>/{kraken2,mapseq}/`. See `training/README.md`.

## 2026-04-30 / 2026-05-01 — V5-V7 amplicon panel

Added the V5-V7 region (799F `AACMGGATTAGATACCCKG` / 1193R `ACGTCATCCCCACCTTCC`,
~395 bp) as a first-class amplicon panel alongside V4, V3-V4, V1-V9 and RRN.
This is the standard plant-microbiome amplicon — the 799F primer mismatches
chloroplast 16S, suppressing host plastid co-amplification.

- `configs/regions.json`: new `V5-V7` block (primers + 300-500 bp length window).
- `scripts/02_generate.py`: `V5-V7` added to `REGIONS_ALL`.
- `src/extract_reads.py`: thin `extract_amplicon()` wrapper around `cutter()` so
  `02_generate.py` can do in-silico PCR for any new region without code edits.
- `scripts/02b_generate_v5v7.py`: helper that mirrors the existing per-replicate
  convention (10 env-seqs files × 5 ART reps × `count=30`) — used to bootstrap
  the V5-V7 mocks from the same per-biome full-length env-seqs files used by V1-V9.
- `data/ref_dbs/V5-V7/`: 330,825 / 337,507 (98%) Greengenes2 backbone sequences
  in-silico-PCR-extracted; per-classifier dirs match the V4/V3-V4/V1-V9
  convention (kraken2/mapseq freshly built; idtaxa/spingo/naive-bayes
  symlinked to the full-length V1-V9 classifiers).
- `data/ref_dbs/V5-V7/taxonova/gg2_v5v7_model.pt`: dedicated 239 MB Taxonova
  model trained with `--smart_full_data`, 99-epoch scout + 99-epoch Stage 2.
- `data/mock_community/<mock>/V5-V7/`: full mock data for all 7 mocks
  (mockrobiota-16S + 6 biomes); ASV counts comparable to V3-V4 (59-94 ASVs).
- `results/<mock>/V5-V7/`: complete benchmark output for all 8 classifiers ×
  7 mocks (56 prediction sets); included in the global
  `results/evaluation/all_metrics.csv` and `by_metric/*.csv`.
- `scripts/v5v7_table1.py`: aggregates species-level Recall/TP/FP per
  (mock × method) into a Markdown sub-table for the manuscript Table 1.

V5-V7 is now identical in every respect to V4/V3-V4 in pipeline behaviour:
the same `02_generate.py` / `03_assign.py` / `04_evaluate.py` flows accept it
as a `--regions` value, the same per-classifier conventions apply, and the
same eval CSVs include it.

## 2026-05-01 — sensitivity + specificity metrics

Added two new evaluation metrics per Reviewer 1 comment 1.9:

- **sensitivity** = taxon-set TP / (TP + FN). Mathematically identical to the
  existing `detection_rate` — both names point to the same calculation but
  `sensitivity` is the standard binary-classification name preferred by the
  reviewer.
- **specificity** = taxon-set TN / (TN + FP), where the universe is the
  unique-taxa count at each level in the per-region reference DB
  (Greengenes2 backbone for V4 / V3-V4 / V5-V7 / V1-V9; GROND for RRN).
  Quantifies how often a classifier correctly refrains from predicting taxa
  absent from the ground truth.

Implementation:
- `scripts/04_evaluate.py`: extended `METRIC_TO_COL`, added a
  `_universe_size()` helper that loads each per-region taxonomy TSV once,
  and computed both metrics inside `_evaluate_combo()`. Module docstring
  expanded with formal definitions.
- `results/evaluation/by_metric/sensitivity.csv` and `specificity.csv` now
  produced alongside the existing per-metric CSVs.
- `results/evaluation/all_metrics.csv` now has 10,880 long-format rows
  (was 8,160) — adds two metric rows per (mock, region, method, param,
  level) combination.
- No benchmark re-run required: metrics are recomputed from existing
  `taxonomy_final.tsv` predictions in a single eval pass.
