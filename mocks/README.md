# Mock communities

Builders for mock community tiers: physical mocks (Mockrobiota-16S, MBARC-26, a gut
mock, a vaginal mock) and biome-simulated mocks (nine biomes with MGnify abundances, and
their rRNA operon counterparts). Each builder writes, per mock and region, the query
sequences and a truth table; `taxobench import-mocks` turns them into a tier
([docs/taxobench.md](../docs/taxobench.md#mock-communities)).

| Path | Content |
|---|---|
| `inputs/` | the mock inputs, gzipped, with `MD5SUMS`: member sequences (`<mock>_seqs.fasta`) and taxonomy (`<mock>_taxa.txt`) |
| `mgnify/` | MGnify abundance retriever and the snapshot behind the biome mocks ([mgnify/README.md](mgnify/README.md)) |
| `scripts/` | build steps `01` (inputs) and `02` (reads, ASVs, truth); helpers in `scripts/lib/` |
| `configs/` | paths (`pipeline_config.json`), primers and lengths (`regions.json`), mock definitions (`mocks.json`) |

## Environment

The builders run in a QIIME 2 environment with ART and Biopython:

```bash
mamba env create -n taxobench-qiime2 --channel-priority flexible \
    --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-linux-conda.yml
mamba install -n taxobench-qiime2 -y -c bioconda -c conda-forge art seqkit biopython joblib
```

## Reference databases

- Greengenes2, 2024.09 release (<https://greengenes2.ucsd.edu/>), for the 16S regions
  V4, V3-V4, V5-V7 and V1-V9.
- GROND, the rRNA operon database built for GTDB release 207, for the operon region (RRN).

Each region has a directory under `data/ref_dbs/` holding its reference sequences and
taxonomy (`sequences.fasta` and `taxonomy.tsv`; `grond_clean.fasta` and `grond_clean.tsv`
for RRN). Classifier indexes may sit beside them, one subdirectory per tool.

## Building

Paths in `configs/pipeline_config.json` resolve against the repository root, so `data/`
sits inside the clone; set `TAXOBENCH_<KEY>` (for example `TAXOBENCH_DATA_DIR`) to put it
elsewhere.

```bash
python mocks/scripts/01a_unpack_mock_inputs.py                     # all mocks into data/
python mocks/scripts/02_generate.py --mock-name mbarc26 --regions V4 V3-V4
taxobench import-mocks --data-dir data --ref-dir data/ref_dbs --out tiers/mocks
```

| Step | Builds |
|---|---|
| `01a_unpack_mock_inputs.py` | places the inputs of `inputs/` in `data/mock_community/input_data/` |
| `01_prepare.py` | collects the Mockrobiota 16S mocks from a Mockrobiota checkout |
| `01b_prepare_extra_mocks.py` | fetches the MBARC-26 16S sequences from NCBI RefSeq (network access) |
| `02_generate.py` | simulates Illumina reads (ART), extracts the region, runs DADA2 and writes the truth |
| `02b_generate_v5v7.py` | the same for V5-V7 |
| `02c_generate_v1v9_direct.py` | full-length V1-V9 cells of the curated mocks |
| `02d_generate_longread_biomes.py` | full-length V1-V9 and RRN cells of the biome mocks |

The truth of each query is the taxonomy of the mock members its sequence matches with
the fewest mismatches (their lowest common ancestor when several tie); an `N` in a
member sequence counts as a mismatch (`--n-wildcard` treats it as a match instead).
The `02c` and `02d` cells are member sequences with random substitutions, so their truth
is known by construction. A new mock needs its two input files in
`data/mock_community/input_data/` and an entry in `configs/mocks.json`; for a new biome,
retrieve its MGnify abundances first with the scripts in `mgnify/`.
