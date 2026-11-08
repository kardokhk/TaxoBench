# Contributing

Issues and pull requests are welcome at <https://github.com/kardokhk/TaxoBench>.
Run the tests before opening a pull request:

```bash
pip install -e ".[test]"
python -m pytest
```

Add a line to `CHANGELOG.md` for any change in behaviour. A change to scoring needs a
known-answer test in `tests/test_known_answer.py`.

## Adding a classifier

1. Write an adapter: a subclass of `taxobench.adapters.base.Adapter` with a `name`, a
   `version`, a `confidence_kind` (`probability`, `score`, `run_share` or `none`) and a
   `parse(raw, **params)` method that returns one `Parsed` record (query ID, lineage
   per rank, confidence per rank) per query. Leave a rank out where the tool makes no
   call. Register it with `@register` in `src/taxobench/adapters/tools.py`, or, from
   another package, under the `taxobench.adapters` entry point.
2. Add a golden test on a raw output of the real tool to
   `tests/test_adapters_golden.py`, with the fixture in `tests/fixtures/adapters/`.
3. Run the tool with a run spec (`docs/taxobench.md`, section Running classifiers).

A tool that already writes a `Feature ID, Taxon, Confidence` table needs no new adapter:
use `taxonomy-tsv`.

## Adding a mock community

1. Write `<mock>_seqs.fasta` and `<mock>_taxa.txt` (`Feature ID`, `Taxon`), by hand
   or with a `01*` script in `mocks/scripts/`. To share the mock, add both files
   gzipped to `mocks/inputs/` and their md5 sums (of the uncompressed files) to
   `mocks/inputs/MD5SUMS`.
2. Register the mock in `mocks/configs/mocks.json` with the regions it is built for.
3. Build the cells (`mocks/README.md`) and convert them with `taxobench import-mocks`.

For a new biome, retrieve its MGnify abundances first with the scripts in `mocks/mgnify/`.

## Adding a tier

Write the five tier tables of the data model (`docs/taxobench.md`, section Data model)
with their sha256 sidecars (`taxobench.contract.write_table` writes both) and check
them with `taxobench validate --tier <dir>`.
