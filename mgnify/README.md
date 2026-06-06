# mgnify — biome-level abundant-taxa retriever

Standalone helper for pulling the most abundant taxa per rank for a given
**MGnify biome**. The output CSVs feed the biome-mock generation step
(picking the top-N genera/species for each environment when constructing
`biome_*` mocks under `data/mock_community/`).

> Not wired into the main pipeline. Run by hand when you want to add a new
> biome or refresh an existing one.

## Files

| File | Purpose |
|---|---|
| `list_biomes.py` | One-shot: scrape the full MGnify biome catalogue and write `biomes.csv`. |
| `fetch_top_taxa.py` | Main script: for one or more biomes, walk samples → runs → analyses → taxonomy and write top-N CSVs per rank. |
| `biomes.csv` | Lineage table (Environment, API_Lineage, Sample_Count, Category). Used by `fetch_top_taxa.py` to resolve biome short-names. |
| `results/<run-tag>/` | Per-run output dir; contains `<biome>__<rank>__topN.csv`. |
| `results/2025-08-20_snapshot/` | Past run covering ~18 biomes (the source for the current `data/mock_community/biome_*` mocks). |

## Usage

```bash
conda activate taxobench-qiime2     # only requires `requests` and `pandas`

# (optional) refresh the lineage table
python list_biomes.py

# fetch a single biome by short-name
python fetch_top_taxa.py --biomes marine --max-samples 100

# multiple biomes at once
python fetch_top_taxa.py --biomes fecal soil oral --max-samples 100 --top-n 100

# a biome that isn't in biomes.csv yet (use the raw lineage)
python fetch_top_taxa.py --lineage 'root:Environmental:Aquatic:Marine:Oceanic' \
                        --label oceanic --max-samples 50
```

Useful flags:

- `--max-samples N` — cap samples fetched per biome (network-bound; ~5–10 s/sample).
- `--top-n N` — top-N taxa per rank to write (default 100).
- `--run-tag NAME` — override the timestamped output dir name.
- `--verbose` — log HTTP back-off / pagination details to stderr.

## Output format

```
results/<run-tag>/<biome>__<rank>__top<N>.csv
```

CSV columns: `rank, taxon, taxa_level, total_abundance` (one row per taxon,
sorted by abundance descending). Ranks emitted: `species, genus, family,
order, class, phylum, kingdom, domain` — but only the ranks present in the
underlying analyses get a file.

## How biome short-names resolve

`fetch_top_taxa.py` reads `biomes.csv` and maps each `Environment` to a key
by lowercasing and replacing spaces / hyphens with underscores. So
`Environment="Human Lung"` becomes the short-name `human_lung`. Run
`python list_biomes.py --top 30` to see the current top biomes.

## API

MGnify v1: `https://www.ebi.ac.uk/metagenomics/api/v1`. The script tries
the SSU taxonomy endpoint first, then falls back to `/taxonomy`,
`/taxonomy/lsu`, `/taxonomy/unite`. 429 responses are retried with
exponential back-off; 404s are silently skipped (some analyses don't have
all four taxonomy endpoints populated).

## Where the output gets used

After fetching, the top-N CSVs are read by the (manual) biome-mock builder
to pick representative species/genera, which then go into
`data/mock_community/biome_<name>/`. That mock is fed to the regular
pipeline (`scripts/02_generate.py` onward) the same way as Mockrobiota.
