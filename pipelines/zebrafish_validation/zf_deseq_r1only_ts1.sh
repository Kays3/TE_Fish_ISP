#!/usr/bin/env bash
# thinkstation1: DESeq2 on the GSE234993 read-1-only counts (DC-16), same
# clutch-blocked design (~ block + role) and frozen env as the main run.
# Waits for NAS r1only/R1ONLY_READY (written by the dgx run, DC-17).
set -uo pipefail
export OPENBLAS_NUM_THREADS=${DESEQ_THREADS:-4} OMP_NUM_THREADS=${DESEQ_THREADS:-4} MKL_NUM_THREADS=${DESEQ_THREADS:-4}
N=/mnt/nas/kaisar/te-oa-full/zf-te-validation-20261001
G=~/workspace/Geneformer_TE
R=~/micromamba/envs/deseq-frozen/bin/Rscript
L=~/workspace/zf_deseq_ts1/r1only; mkdir -p $L/counts
STUDY=ZF_uhrf1_liver_2023
until [ -s $N/r1only/R1ONLY_READY ]; do sleep 300; done
cp $N/r1only/counts/*.cntTable $L/counts/
awk -F'\t' -v s=$STUDY 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; print "sample_id\ttreatment\tblock"; next}
  $c["study_id"]==s {print $c["run_accession"]"\t"$c["role"]"\t"$c["block"]}' \
  $G/config/zebrafish_perturbation_studies.tsv > $L/sheet.tsv
$R $G/scripts/04_bulk_te_deseq.R $L/counts $L/sheet.tsv $L/$STUDY.r1only.tsv > $L/$STUDY.r1only.log 2>&1
rc=$?
mkdir -p $N/r1only/te_de && cp $L/$STUDY.r1only.tsv $L/$STUDY.r1only.log $N/r1only/te_de/ 2>/dev/null
echo "r1only deseq exit=$rc" > $N/r1only/DESEQ_R1ONLY_DONE
