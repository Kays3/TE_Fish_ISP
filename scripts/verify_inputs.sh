#!/usr/bin/env bash
# Verify checksums: repository files and figures always; the large inputs when DATA_ROOT is set.
#   DATA_ROOT=/path/that/contains/te-oa-full scripts/verify_inputs.sh
set -euo pipefail
cd "$(dirname "$0")/.."
. scripts/_sha.sh
bad=0
sha256_check manifests/repo_files.sha256 | grep -v ': OK$' && bad=1 || true
sha256_check manifests/figures.sha256 | grep -v ': OK$' && bad=1 || true
[ $bad = 0 ] && echo "repository files and figures: all OK"
if [ -n "${DATA_ROOT:-}" ]; then
  t="$(mktemp)"
  tail -n +2 manifests/inputs.sha256.tsv | awk -F'\t' -v r="$DATA_ROOT" '{print $1"  "r"/"$3}' > "$t"
  if sha256_check "$t" | grep -v ': OK$'; then bad=1; else echo "large inputs: all OK"; fi
  rm -f "$t"
else
  echo "DATA_ROOT not set: large inputs not checked"
fi
exit $bad
