#!/usr/bin/env bash
# compute host 1, post hoc (DC-18): GSE234993 DESeq2 with unique-mapping % as a
# covariate. MODE=r1only uses NAS r1only/ counts and STAR logs; MODE=pe uses
# the main run's counts and logs (after COUNTS_READY). Two designs each.
set -uo pipefail
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
MODE=${1:?r1only|pe}
N=${DATA_ROOT}/te-oa-full/zf-te-validation-20261001
G=$REPO_ROOT; R=${RSCRIPT:-Rscript}
STUDY=ZF_uhrf1_liver_2023
if [ "$MODE" = r1only ]; then C=$N/r1only/counts; Q=$N/r1only/qc; else C=$N/counts; Q=$N/qc; fi
L=$ZF_DESEQ/mapcov_$MODE; mkdir -p $L/counts
awk -F'\t' -v s=$STUDY 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next} $c["study_id"]==s {print $c["run_accession"]"\t"$c["role"]"\t"$c["block"]}' \
  $G/config/zebrafish_perturbation_studies.tsv > $L/ids.tsv
printf 'sample_id\ttreatment\tblock\tuniq_pct\n' > $L/sheet.tsv
while IFS=$'\t' read -r id role block; do
  cp $C/$id.cntTable $L/counts/
  u=$(grep 'Uniquely mapped reads %' $Q/$id.Log.final.out | cut -f2 | tr -d '%')
  printf '%s\t%s\t%s\t%s\n' "$id" "$role" "$block" "$u" >> $L/sheet.tsv
done < $L/ids.tsv
for d in block noblock; do
  $R $G/scripts/04b_bulk_te_deseq_mapcov.R $L/counts $L/sheet.tsv $L/$STUDY.$MODE.mapcov_$d.tsv $d > $L/$STUDY.$MODE.mapcov_$d.log 2>&1
  echo "$d exit=$?"
done
mkdir -p $N/mapcov && cp $L/*.tsv $L/*.log $N/mapcov/ 2>/dev/null
