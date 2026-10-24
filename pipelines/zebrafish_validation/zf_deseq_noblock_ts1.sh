#!/usr/bin/env bash
# compute host 1: show the clutch-blocking deviation (DC-19). Reruns the three
# clutch-paired studies with ~ role only (block column blanked), same counts,
# env and thread cap as the main run. Output: NAS te_de_noblock/.
set -uo pipefail
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
N=${DATA_ROOT}/te-oa-full/zf-te-validation-20261001
G=$REPO_ROOT; R=${RSCRIPT:-Rscript}
L=$ZF_DESEQ; O=$L/te_de_noblock; mkdir -p $O/sheets
for st in ZF_uhrf1_liver_2020 ZF_uhrf1_liver_2023 ZF_hUHRF1_OE_liver_2023; do
  awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print; next} {$3=""; print}' $L/sheets/$st.tsv > $O/sheets/$st.tsv
  $R $G/scripts/04_bulk_te_deseq.R $L/counts $O/sheets/$st.tsv $O/$st.tsv > $O/$st.log 2>&1
  echo "$st exit=$? $(grep design $O/$st.log)"
done
mkdir -p $N/te_de_noblock && cp $O/*.tsv $O/*.log $N/te_de_noblock/
