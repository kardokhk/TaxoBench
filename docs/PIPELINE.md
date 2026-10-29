# Pipeline reference

The framework is intentionally **modular** — every step is a standalone CLI script
under `scripts/`, backed by importable helpers in `src/`. You can run the entire
matrix end-to-end or just rerun the parts you care about.

## Step-by-step flow

```
01_prepare.py    →  data/mock_community/input_data/<mock>_seqs.fasta + _taxa.txt
                    (skip if already present — mockrobiota-16S, biome_*, gtdb_* are
                     already prepared and symlinked into data/)
        ↓
02_generate.py        (V4 / V3-V4 / V5-V7 ONLY — short Illumina amplicons)
                  →  data/mock_community/<mock>/<region>/source/{paired,single}/
                       rep_seqs_<mock>.{fasta,qza}
                    data/mock_community/<mock>/<region>/trueish-taxonomies.tsv
                    (also already present for everything we ship)

regenerate_env_biome_repseqs.py    (V1-V9 / RRN — long-read regions)
                  →  data/mock_community/<mock>/<region>/source/single/
                       rep_seqs_<mock>.fasta
                    data/mock_community/<mock>/<region>/env-seqs_<mock>{00..09}.fasta
                  see "Long-read regions" below
        ↓
03_assign.py     →  results/<mock>/<region>/<classifier>/<param-id>/taxonomy_final.tsv
        ↓
04_evaluate.py   →  results/evaluation/all_metrics.csv               (tidy long)
                    results/evaluation/by_metric/<metric>.csv        (wide pivots)

Standalone:
05_realworld.py  →  results/real_world/<study>/<classifier>/<param-id>/taxonomy_final.tsv
06_runtime.py    →  results/runtime/<classifier>_runtime.csv
                    results/runtime/all_runtime_results.csv
```

## Mock matrix this framework supports

| Dataset family    | Mock identifier(s)                                                                    | V4 | V3-V4 | V5-V7 | V1-V9 | RRN |
|-------------------|---------------------------------------------------------------------------------------|----|-------|-------|-------|-----|
| Mockrobiota       | `mockrobiota-16S`                                                                     | yes  | yes     | yes     | yes     | —   |
| Biomes (Greengenes2) | `biome_{sputum,soil,skin,lung,digestive_system,mixed}`                             | yes  | yes     | yes     | yes     | —   |
| Biomes (GTDB/GROND)  | `gtdb_{sputum,soil,skin,lung,digestive_system,mixed}`                              | —  | —     | —     | —     | yes   |
| Manufactured mocks | `mbarc26` (Singer 2016, 26 sp)                                                       | yes  | yes     | yes     | yes     | —   |
| Real-world        | `zuo` (V4, 515FB/806RB) · `mantri`, `virtanen` (V3-V4, 341F/785R)                     | yes  | yes     | —     | —     | —   |

V5-V7 (799F `AACMGGATTAGATACCCKG` / 1193R `ACGTCATCCCCACCTTCC`, ~395 bp) is the standard plant-microbiome amplicon — the 799F primer mismatches chloroplast 16S, suppressing host plastid co-amplification. It is wired into the pipeline identically to V4 / V3-V4: same in-silico PCR (`cutter`/`extract_amplicon`), same ART/DADA2 simulation, same per-classifier reference-DB convention. A pre-trained Taxonova V5-V7 model is shipped at `data/ref_dbs/V5-V7/taxonova/gg2_v5v7_model.pt`. IDTAXA / SPINGO / naive-Bayes for V5-V7 follow the same pattern as the other regions and use the full-length V1-V9 classifiers (symlinked).

(RRN is the full rRNA operon and is only available from the GROND/GTDB build, so
the V-region mocks and RRN mocks are different objects — they share biome
identity but use different reference sequences.)

## Long-read regions (V1-V9 + RRN) — different generator

For curated mock-community inputs (mockrobiota-16S, mbarc26),
V1-V9 generation uses `scripts/02c_generate_v1v9_direct.py`, which mutates
the curated species fasta directly (no MGnify, no backbone) — same recipe
parameters as the env-biome regenerator. For env biomes (gg2 backbone /
GROND), see the dedicated long-read regenerator below.

V1-V9 (~1500 bp) and RRN (~5000 bp) **must NOT** go through `02_generate.py`.
ART is an Illumina simulator capped at ~300 bp per read, so any attempt to
generate full-length amplicons through it silently produces 200 bp rep_seqs
that the V1-V9 / RRN models (trained on full-length references) cannot
classify. Per the taxonova manuscript ("mock community generation"):

> "For long-read sequences (V1-V9, full operon), we introduced random
>  mutations at rates of 0.002–0.02 to simulate sequencing errors."

The dedicated entry point is:

```bash
# V1-V9 — uses Greengenes2 backbone (data/ref_dbs/Greengenes2_backbone/)
python /path/to/workdir/taxonova/tx_wip_08/scripts/regenerate_env_biome_repseqs.py \
    --biome biome_marine --region V1-V9
# RRN — uses GROND (data/ref_dbs/RRN/grond_clean.{fasta,tsv})
python /path/to/workdir/taxonova/tx_wip_08/scripts/regenerate_env_biome_repseqs.py \
    --biome gtdb_marine --region RRN
# all 6 env biome cells in one go
python /path/to/workdir/taxonova/tx_wip_08/scripts/regenerate_env_biome_repseqs.py --all
```

Recipe (executed by the script, matching the lung/digestive/skin V1-V9 +
gtdb_lung/digestive RRN files generated 2026-04-26):
1. Read MGnify genus top-100 CSV for the biome
   (`mgnify/results/2025-08-20_snapshot/<biome>_genus_top_100.csv`).
2. For each MGnify genus, find matching gg2-backbone (V1-V9) or GROND (RRN)
   sequences (tolerant of GTDB-style suffixes like `Vibrio_A_672454`).
3. Sample up to 500 sequences weighted by MGnify abundance.
4. For each of 10 replicates, draw a random subsample of ~75 seqs and
   apply per-base mutations at rate uniform in [0.002, 0.02].
5. Write `<biome>/<region>/env-seqs_<biome>{00..09}.fasta` (full length).
6. Dereplicate union by MD5(seq) → write
   `<biome>/<region>/source/single/rep_seqs_<biome>.fasta` (DADA2-style hashed IDs).

Both `02_generate.py` and `slurm/generate_one.sbatch` now error out when
called with `REGION=V1-V9` or `REGION=RRN`, pointing here.

## Step 03: classifier dispatch

`scripts/03_assign.py` dispatches all 8 classifiers in parallel via joblib. The
heavy lifting is in `src/generate_commands.py` (one `build_<classifier>_command()`
per tool). Per-classifier conventions:

| Classifier  | Command                          | Reference inputs                              | Output |
|-------------|----------------------------------|------------------------------------------------|--------|
| taxonova    | `taxonova classification`        | `<region>/taxonova/*model.pt` (+ `.taxonomy.pt`, `.hierarchy.pt`) | `taxonomy_final.tsv` |
| vsearch     | `qiime feature-classifier classify-consensus-vsearch` | `<region>/vsearch/{sequences,taxonomy}.qza`  | via `qiime_helper.py convert_qza_to_tsv` |
| blast       | `qiime feature-classifier classify-consensus-blast`   | `<region>/blast/{blastdb,taxonomy}.qza`      | via `qiime_helper.py` |
| naive-bayes | `qiime feature-classifier classify-sklearn`           | `<region>/naive-bayes/*sklearn*.qza`         | via `qiime_helper.py` |
| idtaxa      | `Rscript run_idtaxa_classification.R`                 | `<region>/idtaxa/*.{rds,RData}`              | `process_idtaxa_output.py` |
| spingo      | `spingo` (in `taxobench-spingo` env)                  | `<region>/spingo/*.fasta`                    | `process_spingo_output.py` (in main env) |
| kraken2     | `kraken2`                                             | `<region>/kraken2/{hash,opts,taxo}.k2d`      | `process_kraken2_output.py` |
| mapseq      | `mapseq`                                              | `<region>/mapseq/reference.{fasta,tax}`      | `process_mapseq_output.py` |

All classifiers default to **single-parameter "default mode"** (one
parameter combo each — see `configs/classifiers_default.json`). To sweep,
duplicate the config with multi-element value lists and pass with
`--params-config <file.json>`.

### Confidence-threshold knobs

| Classifier  | Knob                       | Default |
|-------------|----------------------------|---------|
| taxonova    | `confidence_threshold`     | 0.7     |
| vsearch     | `p-perc-identity`/`p-min-consensus` | 0.80 / 0.51 |
| blast       | same                       | 0.80 / 0.51 |
| naive-bayes | `p-confidence`             | 0.7     |
| idtaxa      | `p-confidence` (0–100)     | 60      |
| spingo      | `bootstrap`                | 10      |
| kraken2     | `confidence`               | 0.0     |
| mapseq      | `cutoff`                   | 0.4     |

## Step 04: evaluation

`scripts/04_evaluate.py` auto-discovers every `<mock>/<region>` under `results/`
that has predictions plus a trueish-taxonomies file, and writes a single
consolidated tree under `results/evaluation/`:

- `all_metrics.csv` — tidy long format with columns
  `mock, region, method, param, level, metric, value`. Load this in pandas/Excel
  and filter however you want.
- `by_metric/<metric>.csv` — one wide pivot per metric (`precision, recall,
  f-score, detection_rate, l1_error, l2_error`), with rows =
  `mock,region,method,param` and columns = `Class, Order, Family, Genus, Species`.

`l1_error` and `l2_error` are abundance-distribution distances between the
truth and prediction at each level. For each level we count occurrences of
each truncated taxon, normalize to relative frequencies, then compute
`L1 = Σ |p_exp(t) − p_obs(t)|` and `L2 = Σ (p_exp(t) − p_obs(t))²` over the
union of taxa (L2 kept as squared, matching the taxobench_11 convention).
Lower is better; 0 = identical distributions.

Levels: `Class, Order, Family, Genus, Species`. The whole tree is rewritten
atomically each run, so re-running after tweaking one classifier's params
just refreshes the affected rows. Predictions under
`results/<mock>/<region>/<classifier>/<param-id>/` are never touched.

No args needed for the full sweep: `python scripts/04_evaluate.py`. Optional
`--mock-name` / `--regions` flags restrict the scan for debugging. No plots.

Other options:

- `--rules {v2,v2-strict,v1,v1-patched}`: scoring rules, `v2` by default. The
  rules are described in README.md ("Scoring rules"). `v1` is the published
  evaluator and, with the v1 truth, reproduces the paper's values.
- `--truth-dir ROOT`: read the trueish files from `ROOT/mock_community/<mock>/<region>/`
  instead of the data directory (for example a v1 or a rebuilt truth tree).
- `--ref-dir ROOT`: reference root holding `<region>/taxonomy.tsv` and
  `RRN/grond_clean.tsv`, which define the universe for specificity. If the
  configured directory does not exist, `data/ref_dbs` of this checkout is used;
  if no universe is found, specificity is written as NaN.
- `--out-dir DIR`: write the evaluation tree to `DIR` instead of
  `results/evaluation/`.

Under `v2` and `v2-strict`, `all_metrics.csv` also holds `jaccard`, `fp_taxa`,
`unassigned_fraction`, `tp`, `fp`, `fn` and, on the species rows, the group
statistics and `strict_species_f`. In those modes L1/L2 and the taxon-set
metrics are computed over named calls only.

## Step 05: real-world

`scripts/05_realworld.py` runs all 8 classifiers on the rep-seq FASTAs in
`data/real_world_input/`: Zuo et al. (`zuo`, V4 / 515FB-806RB, pediatric UC
gut), Mantri et al. (`mantri`, V3-V4 / 341F-785R, soil), and Virtanen et al.
(`virtanen`, V3-V4 / 341F-785R, vaginal). Each study's region is taken from
`configs/mocks.json` (override with `--region`). **No evaluation step** —
there's no ground truth.

Outputs land at `results/real_world/<study>/<classifier>/<param-id>/taxonomy_final.tsv`.

## Step 06: runtime / RAM

`scripts/06_runtime.py` samples a representative fasta at sizes `[100, 500, 1000,
2000, 5000]`, runs each classifier `--replicates` times, and captures wall time
+ peak RSS via `psutil`. Output: one CSV per classifier under `results/runtime/`,
plus `all_runtime_results.csv`.

## Parallelism strategy

Three levels:

1. **Within each classifier**: tools are passed `--threads N` (`--p-threads`,
   `--p-num-threads`, `--p-n-jobs`, `-nthreads`, etc.).
2. **Across classifiers within one mock × region**: `joblib` `Parallel(n_jobs=K)`
   in `03_assign.py`. K=2 is reasonable on 64-core nodes if each classifier uses
   16–20 threads.
3. **Across mocks × regions**: Slurm — submit each (mock × region) as a separate
   sbatch job (see `slurm/assign_full_matrix.sh`).

## Cleaning up to re-run

Anywhere you want a clean re-run, simply remove the destination dir:

```bash
# clean V4 results for one mock then re-run
rm -rf results/mockrobiota-16S/V4
python scripts/03_assign.py --mock-name mockrobiota-16S --regions V4

# same pattern for any region — V4, V3-V4, V5-V7, V1-V9 all behave identically
rm -rf results/mockrobiota-16S/V5-V7
python scripts/03_assign.py --mock-name mockrobiota-16S --regions V5-V7

# clean only one classifier
rm -rf results/mockrobiota-16S/V4/blast
python scripts/03_assign.py --mock-name mockrobiota-16S --regions V4 --classifiers blast
```

Step 02 supports `--force` to wipe previously-generated mock data; step 03
silently overwrites existing per-job output dirs.
