# Data acquisition and layout

None of the large inputs ship with this repo: the reference databases are tens of
gigabytes, the generated mock reads and real-world tables are large, and the
trained models are big binaries. This document tells you what the pipeline expects
on disk, where to get it, and how to lay it out so a run actually goes through.

A note up front: the directory paths in `configs/pipeline_config.json` are
relative and resolve against the repo root, so `data/`, `results/` and so on sit
inside your clone. To put them elsewhere, edit `pipeline_config.json` or set
`TAXOBENCH_<KEY>` (for example `TAXOBENCH_DATA_DIR=/scratch/taxobench/data`).
The conda environment paths still have to be edited by hand. Everything below
is described relative to `data_dir` (the `data/` directory at the repo root by
default). Some scripts also hardcode a `ROOT` path near the top — grep for `ROOT =`
if a region-specific generator complains about a missing path.

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

- **Greengenes2, 2024.09 release** — the backbone for the short/standard 16S
  regions (V4, V3-V4, V5-V7, V1-V9). Download from
  <https://greengenes2.ucsd.edu/>.
- **GROND** (and GTDB-derived rrn sequences) — the full rRNA operon reference for
  the RRN region. See <https://github.com/PaulCotter/GROND>.

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

How the per-tool indexes are located is in `scripts/03_assign.py`
(`_ref_dict`): BLAST/VSEARCH/naive-bayes/IDTAXA each look for their artifact under
`<region>/<tool>/`, Kraken 2 wants `hash.k2d / opts.k2d / taxo.k2d` under
`<region>/kraken2/`, MAPseq wants `reference.fasta` + `reference.tax` under
`<region>/mapseq/`, and SPINGO wants a `.fasta` database. Build these from the
Greengenes2 / GROND sequences with each tool's standard indexing command (QIIME 2
artifacts for the QIIME-based tools, `kraken2-build`, etc.).

### Taxonova models

The per-region Taxonova checkpoints are the `.pt` files above:

| Region | Filename |
|--------|----------|
| V4 | `gg2_v4_mrgroups_model.pt` |
| V3-V4 | `gg2_v3v4_mrgroups_model.pt` |
| V5-V7 | `gg2_v5v7_mrgroups_model.pt` |
| V1-V9 | `gg2_v1v9_mrgroups_model.pt` |
| RRN | `gtdb_rrn_mrgroups_model.pt` |

`scripts/03_assign.py` resolves the Taxonova model by globbing
`<region>/taxonova/*model.pt`, so the filename matters only loosely, but these are
the canonical names. The pretrained checkpoints are published on the Taxonova
repo's Releases page: <https://github.com/kardokhk/taxonova>. Place each region's
`.pt` under `data/ref_dbs/<region>/taxonova/`.

## Tier 1 and 2 mock inputs

Both the physical mocks (tier 1) and the biome-simulated mocks (tier 2) start from
a pair of files per mock in `data/mock_community/input_data/`:

```
<mock>_seqs.fasta     # reference sequences for the mock's members
<mock>_taxa.txt       # their taxonomy
```

The expected filenames per mock are listed in `configs/mocks.json` (`input_seqs`,
`input_taxa`, or `v_seqs`/`v_taxa` for the biomes). For example
`mockrobiota-16S_seqs.fasta`, `mbarc26_seqs.fasta`,
`biome_soil_seqs.fasta`, and so on.

- **Mockrobiota** mocks: `scripts/01_prepare.py` builds these by scanning a
  Mockrobiota checkout at `data/softy/mockrobiota/data/`
  (`mockrobiota_source_dir` in the config). Clone Mockrobiota and point the config
  at it. If you already have prepared `input_data` files symlinked in, step 01 is a
  no-op.
- **MBARC-26**: `scripts/01b_prepare_extra_mocks.py` fetches the 16S records from
  NCBI RefSeq and writes `mbarc26_{seqs,taxa}`.
- **Biome mocks** (`biome_*`): built by sampling Greengenes2 sequences at
  biome-realistic genus abundances. Those abundance profiles come from MGnify and
  are checked in under `mgnify/results/2025-08-20_snapshot/` (top-100 taxa per rank
  per biome). The retriever that produced them is in `mgnify/` — see
  `mgnify/README.md`. Each biome's full-operon counterpart (`gtdb_*`, for RRN) is
  built from GTDB/GROND sequences.

Once the `input_data` pairs exist, `scripts/02_generate.py` (and the `02b`/`02c`
variants) turn them into per-region simulated reads under
`data/mock_community/<mock>/<region>/` and produce the trueish taxonomy that the
evaluator scores against.

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
woltka into an OGU table). `docs/REALWORLD.md` walks through onboarding a new study,
and `docs/REALWORLD_AMPLICON_QC.md` documents the primer/denoising QC. After that,
`scripts/05_realworld_classify.py` and `05_realworld_evaluate.py` are cheap to
rerun.

Study accessions (the three headline studies plus the development studies):

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

1. Set the conda environment paths in `configs/pipeline_config.json` (and the
   directory paths, if `data/` and `results/` should not live in the clone).
2. Put the region's reference DB and per-tool indexes under `data/ref_dbs/<region>/`,
   including the Taxonova `.pt` from the Taxonova Releases page.
3. Get the mock's `input_data` pair in place (clone Mockrobiota, or run the `01*`
   prep script).
4. Run `02_generate.py` -> `03_assign.py` -> `04_evaluate.py` for that
   (mock, region).

For the real-world tier, drop a study's `real_world_input/<study>/` tree in place
and run the two `05_*` scripts.

Be honest with yourself about paths: several scripts were written against the
original absolute layout, so expect to adjust `pipeline_config.json` and the
occasional in-script `ROOT` before a fresh checkout runs clean.
