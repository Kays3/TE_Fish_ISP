#!/usr/bin/env bash
set -euo pipefail
W=~/workspace/zf_te_validation; E=$W/env/te-ocean-acidification; R=$W/reference
cd $R
[ -f genome.fa ] || zcat Danio_rerio.GRCz11.dna.primary_assembly.fa.gz > genome.fa
[ -f genes.gtf ] || zcat Danio_rerio.GRCz11.112.gtf.gz > genes.gtf
[ -f genome.fa.fai ] || $E/bin/samtools faidx genome.fa
[ -f GRCz11_rmsk_TE.gtf ] || $E/bin/python $W/make_te_gtf.py rmsk.txt.gz genome.fa.fai GRCz11_rmsk_TE.gtf
if [ ! -f star_index/SA ]; then
  mkdir -p star_index
  $E/bin/STAR --runMode genomeGenerate --runThreadN 16 --genomeDir star_index \
    --genomeFastaFiles genome.fa --sjdbGTFfile genes.gtf --sjdbOverhang 100 > star_index.log 2>&1
fi
sha256sum genome.fa genes.gtf GRCz11_rmsk_TE.gtf > reference.sha256
echo "reference done" >> $W/logs/STATUS
