#!/usr/bin/env bash
# thinkstation2: run the zebrafish validation analysis once DESeq2 has finished.
set -uo pipefail
N=/mnt/nas/kaisar/te-oa-full/zf-te-validation-20261001
G=~/workspace/Geneformer_TE
until [ -s $N/DESEQ_DONE ]; do sleep 300; done
cd $G && git pull -q --ff-only
mkdir -p runs/zf_validation_20261001/fish05c
cp $G/docs/runs/2026-10-01_fish_bulk_isp/tables/05c_*_per_regulator.tsv runs/zf_validation_20261001/fish05c/
.venv-bulk/bin/python scripts/08_zebrafish_validation.py --te-de $N/te_de --counts $N/counts \
  --samples $N/samples_all.tsv --manifest config/zebrafish_perturbation_studies.tsv \
  --gene-names $N/gene_names.tsv --fish-05c runs/zf_validation_20261001/fish05c \
  --out runs/zf_validation_20261001 > runs/zf_validation_20261001.log 2>&1
echo "08 exit=$?" >> runs/STATS_STATUS
