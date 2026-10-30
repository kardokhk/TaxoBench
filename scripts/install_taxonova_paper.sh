#!/usr/bin/env bash
# Install Taxonova 1.0.0, the version used for the published runs, into the
# active Python environment.
#
# The v1.0.0 tag lists bin/taxonova in setup.py, but the file is not in the tag,
# so `pip install git+...@v1.0.0` fails. This script checks out the tag, adds the
# missing launcher (it calls taxonova.main:main, the same entry point) and installs.
#
# Usage: bash scripts/install_taxonova_paper.sh
#   TAXONOVA_REPO  source repository (default https://github.com/kardokhk/taxonova)
#   TAXONOVA_TAG   tag to install (default v1.0.0)
set -euo pipefail

repo="${TAXONOVA_REPO:-https://github.com/kardokhk/taxonova}"
tag="${TAXONOVA_TAG:-v1.0.0}"
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
python -m pip install .
python -c "import importlib.metadata as m; print('taxonova', m.version('taxonova'))"
