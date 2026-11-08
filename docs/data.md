# Data acquisition and layout

The reference databases, generated mock reads, real-world tables and model
checkpoints are too large for the repository. This document lists what the pipeline
expects on disk, where to obtain it and where to place it.

The directory paths in `configs/pipeline_config.json` are relative and resolve
against the repository root, so `data/`, `results/` and so on sit inside the clone.
To put them elsewhere, edit `pipeline_config.json` or set `TAXOBENCH_<KEY>` (for
example `TAXOBENCH_DATA_DIR=/scratch/taxobench/data`). Set the conda environment
prefixes once, in `conda_envs` of the same file (docs/environments.md). Everything
below is relative to `data_dir` (`data/` at the repository root by default).

The mock community inputs are in this repository (`mocks/inputs/`), with one
published cell (gut mock, V5-V7, classifier outputs and truth) under
`tests/fixtures/` as an example and a regression fixture.

## Top-level data layout

```
data/
  ref_dbs/                     # per-region reference databases + per-tool indexes
  mock_community/
    input_data/                # <mock>_seqs.fasta + <mock>_taxa.txt (tier 1 & 2 inputs)
    <mock>/<region>/...        # generated reads + trueish taxonomy (produced by step 02)
  real_world_input/            # tier 3: rep-seqs, DADA2 tables, shotgun tables per study
  softy/mockrobiota/data/      # Mockrobiota source checkout (tier 1 prep)
```

The exact directory names come from `pipeline_config.json`:
`reference_database_dir`, `mock_input_dir`, `real_world_input_dir`,
`mockrobiota_source_dir`.

## Reference databases

Two source databases back the whole benchmark:

- **Greengenes2, 2024.09 release**: the backbone for the short/standard 16S
  regions (V4, V3-V4, V5-V7, V1-V9). Download from
  <https://greengenes2.ucsd.edu/>.
- **GROND**, the build for GTDB release 207: the full rRNA operon reference for the
  RRN region (<https://zenodo.org/records/10889037>).

Each region gets its own subdirectory under `data/ref_dbs/`, and within it a common
set of files plus one subfolder per classifier that needs a prebuilt index. Observed
layout:

```
data/ref_dbs/
  V4/
    sequences.fasta           # region reference sequences
    taxonomy.tsv              # reference taxonomy
    merged.fasta             # (used by some tools)
    sequences.qza, taxonomy.qza
    blast/  vsearch/  naive-bayes/  idtaxa/  kraken2/  mapseq/  spingo/
    taxonova/
      gg2_v4_mrgroups_model.pt
  V3-V4/   ... taxonova/gg2_v3v4_mrgroups_model.pt
  V5-V7/   ... taxonova/gg2_v5v7_mrgroups_model.pt
  V1-V9/   ... taxonova/gg2_v1v9_mrgroups_model.pt
  RRN/
    grond_clean.fasta
    grond_clean.tsv
    ... taxonova/gtdb_rrn_mrgroups_model.pt
```

How the per-tool indexes are located is in `scripts/_common.py`
(`ref_dict`): BLAST/VSEARCH/naive-bayes/IDTAXA each look for their artifact under
`<region>/<tool>/`, Kraken 2 wants `hash.k2d / opts.k2d / taxo.k2d` under
`<region>/kraken2/`, MAPseq wants `reference.fasta` + `reference.tax` under
`<region>/mapseq/`, and SPINGO wants a `.fasta` database. Build these from the
Greengenes2 / GROND sequences with each tool's standard indexing command (QIIME 2
artifacts for the QIIME-based tools, `kraken2-build`, etc.).

### Taxonova models

The default panel (`configs/classifiers_default.json`) uses model ids, and the
Taxonova package downloads and verifies the models itself:

```bash
taxonova models download all        # into $TAXONOVA_MODELS_DIR, default ~/.taxonova/models
```

The paper configuration (`configs/classifiers_paper.json`) uses the v1.0.0 checkpoints
from the Taxonova v1.0.0 release, placed under `data/ref_dbs/<region>/taxonova/`:

| Region | File |
|--------|------|
| V4 | `gg2_v4_mrgroups_model.pt` |
| V3-V4 | `gg2_v3v4_mrgroups_model.pt` |
| V5-V7 | `gg2_v5v7_mrgroups_model.pt` |
| V1-V9 | `gg2_v1v9_mrgroups_model.pt` |
| RRN | `gtdb_rrn_mrgroups_model.pt` |

Each entry of that config gives the download URL and md5, and `03_assign.py` checks the
md5 before a run. A config without a `models` entry for a region uses the first
`<region>/taxonova/*model.pt`.

## Tier 1 and 2 mock inputs

Both the physical mocks (tier 1) and the biome-simulated mocks (tier 2) start from
a pair of files per mock in `data/mock_community/input_data/`:

```
<mock>_seqs.fasta     # reference sequences for the mock's members
<mock>_taxa.txt       # their taxonomy (Feature ID, Taxon)
```

The files used for the published cells are in this repository, gzipped, in
`mocks/inputs/` (with `MD5SUMS`): Mockrobiota-16S, MBARC-26, the gut and vaginal mocks,
the nine `biome_*` mocks and their nine `gtdb_*` RRN counterparts. Place them with

```bash
python scripts/01a_unpack_mock_inputs.py              # all mocks
python scripts/01a_unpack_mock_inputs.py --mock-name mbarc26 biome_soil
```

The expected file names per mock are listed in `configs/mocks.json` (`input_seqs`,
`input_taxa`, or `v_seqs`/`v_taxa` for the biomes). How the inputs were built:

- **Mockrobiota-16S**: the 16S mocks of a Mockrobiota checkout
  (`mockrobiota_source_dir`), collected by `scripts/01_prepare.py`.
- **MBARC-26**: 16S records from NCBI RefSeq, fetched by
  `scripts/01b_prepare_extra_mocks.py` (needs Biopython and network access). The
  published taxonomy labels were then resolved against the Greengenes2 2024.09
  taxonomy, so they differ from the labels 01b writes; use the file in
  `mocks/inputs/` to reproduce the published cells.
- **Biome mocks** (`biome_*`): Greengenes2 sequences sampled at biome-specific genus
  abundances from MGnify. The abundance profiles are in
  `mocks/mgnify/results/2025-08-20_snapshot/` (top-100 taxa per rank per biome), retrieved
  with the scripts in `mocks/mgnify/` (see `mocks/mgnify/README.md`). Each biome's full-operon
  counterpart (`gtdb_*`, for RRN) was built from GROND sequences. The sampling
  script is not part of this repository; the inputs it produced are.

`scripts/02_generate.py` (and the `02b`, `02c` and `02d` variants) turn the inputs
into per-region cells under `data/mock_community/<mock>/<region>/` and write the
truth that the evaluator scores against.

## Tier 3 real-world studies

Each paired study lives under `data/real_world_input/<study>/` with the layout
referenced in `configs/mocks.json` (`real_world` block):

```
real_world_input/
  rep_seqs_<study>.fasta            # amplicon ASV sequences (classifier input)
  <study>/
    dada2_output/asv_table_ids.csv  # amplicon ASV count table
    woltka_output/ogu_table.biom    # paired shotgun OGU (operational genome unit) table
    qiime2_output/exported_taxonomy/taxonomy.tsv   # OGU -> taxonomy
    metadata/paired_samples.csv     # amplicon_16s_run <-> shotgun_wgs_run pairing
```

You generate these once per study by downloading the raw reads, primer-trimming and
denoising the amplicon side (DADA2), and profiling the shotgun side (bowtie2 +
woltka into an OGU table). `docs/realworld.md` walks through onboarding a new study,
and `docs/realworld_amplicon_qc.md` documents the primer/denoising QC. After that,
`scripts/05_realworld_classify.py` and `05_realworld_evaluate.py` are cheap to
rerun.

Study accessions. The paper reports the first three (Zuo, Mantri, Virtanen); the
other three were used during development:

| Study | Biome | Region | Amplicon accession | Shotgun accession |
|-------|-------|--------|--------------------|-------------------|
| Zuo et al. | gut (pediatric UC) | V4 | NCBI PRJNA759642 | (same project) |
| Mantri et al. | soil | V3-V4 | NCBI PRJNA717918 | NCBI PRJNA717813 |
| Virtanen et al. | vaginal | V3-V4 | ENA PRJEB25778 | ENA PRJEB43447 |
| Cao et al. (coral) | Florida coral reef | V4 | NCBI PRJNA910824 | (same project) |
| EMP500 forest soil (Myrold) | soil | V4 | PRJEB42019 (sub-study myrold.59) | PRJEB42019 |
| EMP500 plant phyllosphere (Berry) | plant | V4 | PRJEB42019 (sub-study berry.2) | PRJEB42019 |

For Zuo, the paired-sample mapping was assembled from study metadata (see the note
in `configs/mocks.json`); the EMP500 entries are specific sub-studies of the larger
EMP500 project (PRJEB42019).

## Putting it together

A minimal reproduction for one region, one mock would be:

1. Set the conda environment prefixes in `configs/pipeline_config.json` (and the
   directory paths, if `data/` and `results/` should not live in the clone).
2. Put the region's reference DB and per-tool indexes under `data/ref_dbs/<region>/`,
   and download the Taxonova models (`taxonova models download all`).
3. Place the mock's input pair (`scripts/01a_unpack_mock_inputs.py`).
4. Run `02_generate.py` -> `03_assign.py` -> `04_evaluate.py` for that
   (mock, region).

For the real-world tier, drop a study's `real_world_input/<study>/` tree in place
and run the two `05_*` scripts.
