# Contributing

Issues and pull requests are welcome at <https://github.com/kardokhk/TaxoBench>.
Run the tests before opening a pull request:

```bash
cd tests && python -m unittest discover -s . -p "test_*.py"
```

Add a line to `CHANGELOG.md` for any change in behavior. Changes to scoring must keep
`--rules v1` reproducing the published values (`tests/test_regression_v1.py`, and
`tests/test_golden_v1.py` over all published cells).

## Adding a classifier

1. Add a parser `src/process_<tool>_output.py` that converts the tool's output to
   `taxonomy_final.tsv` with the columns `Feature ID`, `Taxon` and `Confidence`.
   `Taxon` is a lineage `d__...; p__...; c__...; o__...; f__...; g__...; s__...`;
   leave a rank empty (`s__`) when the tool makes no call at that rank.
2. Add `build_<tool>_command()` to `src/generate_commands.py`, returning one shell
   command that writes that file, and register it in the `builders` table and the
   `SUPPORTED` tuple of `scripts/03_assign.py` (and `05_realworld_classify.py` for the
   real-world tier).
3. If the tool needs a prebuilt index, add its lookup to `ref_dict()` in
   `scripts/_common.py` (indexes live under `data/ref_dbs/<region>/<tool>/`).
4. Add the tool to `classifiers`, `method_conda_envs` and
   `method_parameters_combinations` in `configs/classifiers_default.json`.
5. Add a test in `tests/test_scoring_rules.py` for the tool's abstention format if it
   differs from those already covered.

To score outputs produced elsewhere, none of these steps is needed: place the files as
`<results>/<mock>/<region>/<tool>/<parameters>/taxonomy_final.tsv` and run
`scripts/04_evaluate.py --results-dir <results>`.

## Adding a mock community

1. Write `<mock>_seqs.fasta` and `<mock>_taxa.txt` (`Feature ID`, `Taxon`), by hand
   or with a `01*` script. To share the mock, add both files gzipped to
   `mocks/inputs/` and their md5 sums (of the uncompressed files) to
   `mocks/inputs/MD5SUMS`; `01a_unpack_mock_inputs.py` then places them in the input
   directory (`mock_input_dir` in `configs/pipeline_config.json`).
2. Register the mock in `configs/mocks.json` with the regions it is evaluated on.
3. Generate the cells: `02_generate.py` for V4, V3-V4 and V5-V7;
   `02c_generate_v1v9_direct.py` for V1-V9; `02d_generate_longread_biomes.py` for
   long-read cells of biome mocks.
4. Run `03_assign.py` and `04_evaluate.py` for the new mock.

For a new biome, retrieve its MGnify abundances first with the scripts in `mgnify/`.

## Adding a Taxonova model

New Taxonova models (other references or regions) are published by the Taxonova
package under a model id. Add the id for the region to `models` in the `taxonova` block
of `configs/classifiers_default.json`, or create a separate classifier config that
points to it.
