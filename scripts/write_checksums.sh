#!/usr/bin/env bash
# Rewrite manifests/repo_files.sha256 and manifests/figures.sha256 from the working tree.
set -euo pipefail
cd "$(dirname "$0")/.."
. scripts/_sha.sh
LC_ALL=C find README.md requirements.txt data config docs scripts pipelines environment -type f -not -path '*/__pycache__/*' | sort | while read -r f; do sha256_make "$f"; done > manifests/repo_files.sha256
LC_ALL=C find docs/runs -path '*/figures/*.png' -type f | sort | while read -r f; do sha256_make "$f"; done > manifests/figures.sha256
echo "wrote $(wc -l < manifests/repo_files.sha256) + $(wc -l < manifests/figures.sha256) checksums"
