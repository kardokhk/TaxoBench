#!/usr/bin/env bash
# Install Taxonova 1.0.0, the version used for the published runs, with the
# dependency versions of the environment that ran them, into the active Python
# environment (Python 3.10 recommended; the paper runs used 3.10.14).
#
# Two defects of the v1.0.0 tag are worked around here:
#  - setup.py lists bin/taxonova, but the file is not in the tag, so
#    `pip install git+...@v1.0.0` fails. The launcher is added (it calls
#    taxonova.main:main, the same entry point).
#  - pyproject.toml has a [project] table without dependencies, which overrides
#    setup.py's install_requires, so `pip install .` installs Taxonova with no
#    dependencies at all and `taxonova` then fails on `import setproctitle`.
#    The dependencies are therefore installed explicitly, pinned to the
#    versions recorded for the paper environment (taxobench-qiime2, 2026-09-14),
#    and the package itself with --no-deps.
#
# The script ends with a check that the command starts; it exits non-zero on
# any failure.
#
# Usage: bash scripts/install_taxonova_paper.sh
#   TAXONOVA_REPO         source repository (default https://github.com/kardokhk/taxonova)
#   TAXONOVA_TAG          tag to install (default v1.0.0)
#   TAXONOVA_TORCH_INDEX  wheel index for torch (default the CPU index
#                         https://download.pytorch.org/whl/cpu; the paper runs
#                         used the CUDA 12.8 build, cu128, and inference on CPU)
set -euo pipefail

repo="${TAXONOVA_REPO:-https://github.com/kardokhk/taxonova}"
tag="${TAXONOVA_TAG:-v1.0.0}"
torch_index="${TAXONOVA_TORCH_INDEX:-https://download.pytorch.org/whl/cpu}"

# Versions of the environment that produced the published outputs
# (Python 3.10.14; torch was pip-installed, the rest came from conda-forge).
pins=(
    "numpy==1.26.4"
    "pandas==2.2.2"
    "python-dateutil==2.8.2"
    "pytz==2024.2"
    "scikit-learn==1.4.2"
    "scipy==1.13.0"
    "joblib==1.4.2"
    "threadpoolctl==3.5.0"
    "tqdm==4.66.5"
    "psutil==6.1.0"
    "setproctitle==1.3.3"
)

py_version="$(python -c 'import sys; print("%d.%d.%d" % sys.version_info[:3])')"
case "$py_version" in
    3.10.*) ;;
    *) echo "warning: Python $py_version; the paper runs used 3.10.14 and the pins above are for 3.10" >&2 ;;
esac

work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

git clone --quiet --depth 1 --branch "$tag" "$repo" "$work/taxonova"
cd "$work/taxonova"
if [ ! -f bin/taxonova ]; then
    mkdir -p bin
    printf '%s\n' '#!/usr/bin/env python' \
        'from taxonova.main import main' \
        '' \
        'if __name__ == "__main__":' \
        '    main()' > bin/taxonova
    chmod +x bin/taxonova
fi

python -m pip install "torch==2.7.1" --index-url "$torch_index"
python -m pip install "${pins[@]}"
python -m pip install --no-deps .

python -c "import importlib.metadata as m; print('taxonova', m.version('taxonova'), '| torch', m.version('torch'), '| numpy', m.version('numpy'), '| pandas', m.version('pandas'))"
# Fails (and with set -e stops the script) if any import is missing.
taxonova classification --help > /dev/null
echo "taxonova classification starts"
