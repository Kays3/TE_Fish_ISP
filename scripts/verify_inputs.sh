#!/usr/bin/env bash
# Verify checksums: repository files and figures always; the large inputs when DATA_ROOT is set.
#   DATA_ROOT=/path/that/contains/te-oa-full scripts/verify_inputs.sh
# Exit status 1 if any file is missing or differs.
set -uo pipefail
cd "$(dirname "$0")/.."
. scripts/_sha.sh
bad=0
check() {  # check LABEL CHECKFILE
  local out
  if out="$(sha256_check "$2" 2>&1)"; then echo "$1: all $(wc -l < "$2" | tr -d ' ') OK"; else
    echo "$1: PROBLEM"; echo "$out" | grep -v ': OK$' | head -n 10; bad=1; fi
}
check "repository files" manifests/repo_files.sha256
check "committed figures" manifests/figures.sha256
if [ -n "${DATA_ROOT:-}" ]; then
  t="$(mktemp)"
  tail -n +2 manifests/inputs.sha256.tsv | awk -F'\t' -v r="$DATA_ROOT" '{print $1"  "r"/"$3}' > "$t"
  check "large inputs" "$t"
  rm -f "$t"
else
  echo "DATA_ROOT not set: large inputs not checked"
fi
exit $bad
