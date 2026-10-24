#!/usr/bin/env bash
# One-command driver for the fish ISP work.
#
#   scripts/run_all.sh            # = figures: rebuild the 10 figures from docs/runs/*/tables and compare with the committed PNGs (seconds, no data)
#   scripts/run_all.sh verify     # checksums of repository files, figures and (DATA_ROOT set) the large inputs
#   scripts/run_all.sh compute    # recompute the tables from the counts (needs DATA_ROOT; heavy, see below)
#   scripts/run_all.sh all        # verify, compute, figures
#
# Environment: PYTHON (default python3; use the pinned environment, environment/README.md),
#              DATA_ROOT (directory that contains te-oa-full/), THREADS (default 10).
#
# `compute` reruns the CPU analyses 05 to 08 on the archived counts. It is NOT part of the first
# migration check: the commands are copied from the run records of 1 to 2 October 2026 and were not
# re-executed when this repository was made. Step 05 took about 5 minutes (primary) to 20 minutes
# (axis removed) on 20 cores; the random-regulator runs took about 12 minutes each on 10 threads.
# The zebrafish alignment and DESeq2 steps (pipelines/zebrafish_validation/) are not rerun here: their
# outputs are the manifested inputs counts/, te_de/ and qc/. No GPU is used anywhere.
set -euo pipefail
cd "$(dirname "$0")/.."
PY="${PYTHON:-python3}"
THREADS="${THREADS:-10}"

stage_verify() { scripts/verify_inputs.sh; }

stage_figures() { "$PY" scripts/make_figures.py --out rebuilt --check; }

stage_compute() {
  : "${DATA_ROOT:?set DATA_ROOT to the directory that contains te-oa-full/}"
  export DATA_ROOT OMP_NUM_THREADS="$THREADS" OPENBLAS_NUM_THREADS="$THREADS"
  mkdir -p runs
  cp data/small_inputs/*.tsv runs/
  # configs: the historical NAS prefix (/mnt/nas/kaisar) is replaced by DATA_ROOT; runs/ is the working directory the configs refer to
  for f in docs/runs/2026-10-01_fish_bulk_isp/configs/*.yaml docs/runs/2026-10-01_gene_axis/configs/*.yaml \
           docs/runs/2026-10-01_te_axis_treatment/configs/*.yaml; do
    sed "s|/mnt/nas/kaisar|$DATA_ROOT|g" "$f" > "runs/$(basename "$f")"
  done
  for s in 2 3; do
    sed "s/fish_random1/fish_random$s/g" runs/fish_random1.yaml > "runs/fish_random$s.yaml"
    sed "s/fish_random1/fish_random$s/g" runs/fish_random1_rmTEpc1.yaml > "runs/fish_random${s}_rmTEpc1.yaml"
  done

  # 05: bulk linear ISP (CellOracle-style), real regulators, then expression-matched random regulators
  for r in fish_bulk_isp_20261001 fish_bulk_isp_20261001_no_outliers fish_bulk_isp_20261001_rmTEpc1; do
    "$PY" scripts/05_bulk_linear_isp.py --config "runs/$r.yaml"
  done
  for s in 1 2 3; do
    "$PY" scripts/05_bulk_linear_isp.py --config "runs/fish_random$s.yaml" --random-regulators "$s"
    "$PY" scripts/05_bulk_linear_isp.py --config "runs/fish_random${s}_rmTEpc1.yaml" --random-regulators "$s"
  done
  # 05b: the global TE-load axis; 05c: real regulators against random pseudo-regulators
  "$PY" scripts/05b_te_axis_diagnostics.py --config runs/fish_bulk_isp_20261001.yaml --qc runs/qc.tsv
  R=runs/fish_bulk_isp_20261001/tables; Rr=runs/fish_bulk_isp_20261001_rmTEpc1/tables
  "$PY" scripts/05c_compare_random_controls.py --real "$R/05_bulk_isp_effects.tsv" \
      --random runs/fish_random{1,2,3}/tables/05_random{1,2,3}_bulk_isp_effects.tsv --out runs/05c --label plain
  "$PY" scripts/05c_compare_random_controls.py --real "$Rr/05_rmTEpc1_bulk_isp_effects.tsv" \
      --random runs/fish_random{1,2,3}_rmTEpc1/tables/05_random{1,2,3}_rmTEpc1_bulk_isp_effects.tsv --out runs/05c --label rmTEpc1
  # 06 and 07: does CO2 shift the TE axis, and which genes account for it
  "$PY" scripts/06_te_axis_treatment.py --config runs/fish_te_axis_treatment_20261001.yaml \
      --qc runs/qc_align.tsv --de runs/oa_significant_te_results.tsv
  "$PY" scripts/07_treatment_gene_axis.py --config runs/fish_gene_axis_20261001.yaml \
      --gene-sets config/gene_sets_apoly.yaml --n-random-sets 10000 --n-boot 1000 --n-random-genes 3000
  # 08: zebrafish ground truth, from the manifested DESeq2 tables, counts and QC
  N="$DATA_ROOT/te-oa-full/zf-te-validation-20261001"
  mkdir -p runs/zf_validation_20261001/fish05c
  cp runs/05c/05c_*_per_regulator.tsv runs/zf_validation_20261001/fish05c/
  "$PY" scripts/08_zebrafish_validation.py --te-de "$N/te_de" --counts "$N/counts" \
      --samples "$N/samples_all.tsv" --manifest config/zebrafish_perturbation_studies.tsv \
      --gene-names "$N/gene_names.tsv" --fish-05c runs/zf_validation_20261001/fish05c --out runs/zf_validation_20261001
  T=runs/zf_validation_20261001/tables
  "$PY" scripts/08c_qc_by_arm.py "$N/qc" "$N/samples_all.tsv" "$T"
  "$PY" scripts/08d_te_share_contributors.py "$N/counts" "$N/samples_all.tsv" ZF_hUHRF1_OE_liver_2023 \
      "$T/08d_hUHRF1_OE_te_share_contributors.tsv" "$N/qc"
  "$PY" scripts/08f_exploratory_concordance.py "$T"
  echo "compute done; new tables are under runs/ (compare with docs/runs/*/tables using diff or the 05c/08 tables)"
  echo "08b and 08e need the DESeq2 working folders (mapcov, r1only, noblock): see pipelines/zebrafish_validation/"
}

case "${1:-figures}" in
  figures) stage_figures ;;
  verify) stage_verify ;;
  compute) stage_compute ;;
  all) stage_verify; stage_compute; stage_figures ;;
  *) echo "usage: $0 [figures|verify|compute|all]"; exit 2 ;;
esac
