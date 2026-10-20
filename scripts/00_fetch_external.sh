#!/usr/bin/env bash
# Fetch reproducible external data. TCRvdb is request-only and NOT fetched here.
set -euo pipefail
cd "$(dirname "$0")/../data/external"
VER=2026-06-03
if [ ! -d "vdjdb/vdjdb-${VER}" ]; then
  curl -fsSL -o "vdjdb-${VER}.zip" \
    "https://github.com/antigenomics/vdjdb-db/releases/download/${VER}-ZENODO/vdjdb-${VER}.zip"
  unzip -oq "vdjdb-${VER}.zip" -d vdjdb
fi
sha256sum "vdjdb/vdjdb-${VER}/vdjdb_full.txt" | tee vdjdb_full.sha256
echo "TCRvdb: request at https://tcrvdb.nkiai.nl -> place CSV in data/external/tcrvdb/"
