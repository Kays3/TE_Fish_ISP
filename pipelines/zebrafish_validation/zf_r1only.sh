#!/usr/bin/env bash
# pSS: sensitivity analysis for GSE234993 (ZF_uhrf1_liver_2023). Mutants map
# ~7 points lower than controls (unique %), driven by "too short" pairs after
# adapter trimming, so realign that study with read 1 only. Waits for the main
# alignment to exit, because align-and-count.sh removes STAR's shared genome on
# exit. Outputs to rnaseq_r1only/ and NAS r1only/.
set -uo pipefail
W=~/workspace/zf_te_validation; E=$W/env/te-ocean-acidification; R=$W/reference
N=/mnt/nas/kaisar/te-oa-full/zf-te-validation-20261001
export PATH=$W/wrappers:$E/bin:$PATH
until grep -q "align exit" $W/logs/STATUS 2>/dev/null; do sleep 300; done
STUDY=ZF_uhrf1_liver_2023
awk -F'\t' -v s=$STUDY 'NR==1 || $5==s {OFS="\t"; if (NR>1) $3="-"; print}' $W/samples_all.tsv > $W/samples_r1only.tsv
THREADS=8 PARALLEL_JOBS=3 TECOUNT_JOBS=9 BAM_SORT_RAM=30000000000 \
  bash $W/align-and-count.sh $W/samples_r1only.tsv $R/star_index $R/genes.gtf $R/GRCz11_rmsk_TE.gtf $W/rnaseq_r1only >> $W/logs/r1only.log 2>&1
rc=$?
mkdir -p $N/r1only/counts $N/r1only/qc
for f in $W/rnaseq_r1only/counts/*/*.cntTable; do cp "$f" $N/r1only/counts/; done
for f in $W/rnaseq_r1only/bam/*.Log.final.out; do cp "$f" $N/r1only/qc/; done
n=$(find $N/r1only/counts -name '*.cntTable' -size +0c | wc -l)
echo "r1only exit=$rc nas=$n" >> $W/logs/STATUS
echo "r1only exit=$rc nas=$n" > $N/r1only/R1ONLY_READY
