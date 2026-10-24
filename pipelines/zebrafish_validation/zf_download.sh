#!/usr/bin/env bash
# Resumable, checksum-verified download of the zebrafish validation panel (ENA + references).
set -uo pipefail
W=$ZF_WORK
mkdir -p $W/fastq $W/reference $W/logs
cd $W
fetch() {  # url md5 bytes
  local url=$1 md5=$2 out=$W/fastq/$(basename "$1")
  [ -f "$out.ok" ] && return 0
  for attempt in 1 2 3 4 5; do
    curl -sS -L --retry 5 --retry-delay 20 -C - -o "$out" "$url" && \
    [ "$(md5sum "$out" | cut -d' ' -f1)" = "$md5" ] && { touch "$out.ok"; echo "OK $out"; return 0; }
    echo "retry $attempt $out" >&2; [ $attempt -ge 3 ] && rm -f "$out"; sleep 30
  done
  echo "FAILED $url" >&2; return 1
}
export -f fetch; export W
# references first
cd $W/reference
for u in https://ftp.ensembl.org/pub/release-112/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz \
         https://ftp.ensembl.org/pub/release-112/gtf/danio_rerio/Danio_rerio.GRCz11.112.gtf.gz \
         https://hgdownload.soe.ucsc.edu/goldenPath/danRer11/database/rmsk.txt.gz; do
  [ -f "$(basename $u).ok" ] || { curl -sS -L --retry 5 -C - -O "$u" && touch "$(basename $u).ok"; }
done
curl -sS -L https://ftp.ensembl.org/pub/release-112/fasta/danio_rerio/dna/CHECKSUMS -o CHECKSUMS.ensembl.fasta
curl -sS -L https://ftp.ensembl.org/pub/release-112/gtf/danio_rerio/CHECKSUMS -o CHECKSUMS.ensembl.gtf
cd $W
awk -F'\t' '{print $1" "$2" "$3}' zf_download_list.tsv | xargs -P 6 -n 3 bash -c 'fetch "$0" "$1" "$2"' > logs/download.log 2> logs/download.err
echo "download exit=$? ok=$(ls fastq/*.ok | wc -l) of $(wc -l < zf_download_list.tsv)" >> logs/STATUS
