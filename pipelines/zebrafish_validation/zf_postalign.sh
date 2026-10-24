#!/usr/bin/env bash
# pSS: after alignment, publish the small count outputs to the NAS for DESeq2 on thinkstation1.
set -uo pipefail
W=~/workspace/zf_te_validation
N=/mnt/nas/kaisar/te-oa-full/zf-te-validation-20261001
until grep -q "align exit" $W/logs/STATUS 2>/dev/null; do sleep 300; done
mkdir -p $N/counts $N/qc
for f in $W/rnaseq/counts/*/*.cntTable; do cp "$f" $N/counts/; done
for f in $W/rnaseq/bam/*.Log.final.out; do cp "$f" $N/qc/; done
cp $W/rnaseq/qc/*.fastp.json $N/qc/ 2>/dev/null
zcat $W/reference/Danio_rerio.GRCz11.112.gtf.gz | awk -F'\t' '$3=="gene"' \
  | grep 'gene_name' | sed -E 's/.*gene_id "([^"]+)".*gene_name "([^"]+)".*/\1\t\2/' > $N/gene_names.tsv
cp $W/samples_all.tsv $W/zebrafish_perturbation_studies.tsv $W/reference/reference.sha256 $N/
n_local=$(ls $W/rnaseq/counts/*/*.cntTable | wc -l); n_nas=$(find $N/counts -name '*.cntTable' -size +0c | wc -l)
echo "local=$n_local nas=$n_nas" > $N/COUNTS_READY
echo "postalign done local=$n_local nas=$n_nas" >> $W/logs/STATUS
