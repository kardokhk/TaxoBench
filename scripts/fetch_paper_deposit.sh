#!/usr/bin/env bash
# Download the bundles of the paper's Zenodo deposit needed to rescore the
# published tables (doi:10.5281/zenodo.23808921, version 1.0), check their
# SHA-256 against the deposit's manifest values below, and unpack them.
#
# Usage: bash scripts/fetch_paper_deposit.sh DEST
# Then:  TAXOBENCH_PAPER_DEPOSIT=DEST python -m unittest tests/test_golden_v1.py
set -euo pipefail

dest="${1:?usage: fetch_paper_deposit.sh DEST}"
record="https://zenodo.org/records/23808921/files"
mkdir -p "$dest/_download"

# sha256 values from MANIFEST.sha256 of the deposit
while read -r sha name; do
    fp="$dest/_download/$name"
    if [ ! -f "$fp" ] || ! echo "$sha  $fp" | sha256sum --check --status; then
        curl -fsSL --retry 3 -o "$fp" "$record/$name?download=1"
    fi
    echo "$sha  $fp" | sha256sum --check
    tar -xzf "$fp" -C "$dest"
done <<'EOF'
b3c3aa834282509224e860ab502db2b7c63df7909e63cb41ee41b973478de15e paper_classifier_outputs.tar.gz
c5f18092d55b08fc8450178e370fb23d973fa37964efd00d57132787db8955ed paper_metric_tables.tar.gz
7b5e0248506b823180b3d8bfeeffa0048946dca7fdca99a45ee2b2b924439a4b paper_truth_v1.tar.gz
aa943d5dfd4f0ef338a10b5b5b09cacc34fbc65fa3d6e6783278d8f5b1d56cf3 reference_taxonomies.tar.gz
EOF
echo "unpacked into $dest"
