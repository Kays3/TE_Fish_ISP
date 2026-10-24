#!/usr/bin/env bash
# One-shot status of the zebrafish validation chain, run from a client that
# reaches the workers through the gateway (see git-sync-bastion-trace).
G=${GATEWAY:-gateway-hokudaibook}
# GATEWAY=none connects to the workers directly (works where ~/.ssh/config has them).
if [ "$G" = none ]; then r() { ssh -n -o BatchMode=yes -o ConnectTimeout=15 "$1" "$2" 2>&1; }
else r() { ssh -n -o BatchMode=yes -o ConnectTimeout=15 "$G" "ssh -n -o BatchMode=yes -o ConnectTimeout=15 $1 '$2'" 2>&1; }; fi
echo "== pSS: download / reference / align"
r pSS 'W=~/workspace/zf_te_validation; date -u; cat $W/logs/STATUS 2>/dev/null; echo "fastq verified: $(ls $W/fastq/*.ok 2>/dev/null | wc -l)/220  ($(du -sh $W/fastq | cut -f1))"; echo "failed downloads: $(grep -c FAILED $W/logs/download.err 2>/dev/null)"; echo "counted samples: $(ls $W/rnaseq/counts/*.tecount.complete 2>/dev/null | wc -l)/112  failed: $(ls $W/rnaseq/counts/*.tecount.failed 2>/dev/null | wc -l)"; echo "watchers alive: $(pgrep -f "[z]f_(download|align|postalign)" | wc -l)"; df -h ~ | tail -1; grep r1only $W/logs/STATUS; echo "r1only watcher alive: $(pgrep -f "[z]f_r1only" | wc -l)"'
echo "== NAS markers"
r thinkstation2 'ls /mnt/nas/kaisar/te-oa-full/zf-te-validation-20261001/ 2>/dev/null; cat /mnt/nas/kaisar/te-oa-full/zf-te-validation-20261001/{COUNTS_READY,DESEQ_DONE} 2>/dev/null'
echo "== thinkstation1: DESeq2"
r thinkstation1 'echo "watcher alive: $(pgrep -f "[z]f_deseq_ts1" | wc -l)"; ls ~/workspace/zf_deseq_ts1/te_de/*.tsv 2>/dev/null | wc -l; cat ~/workspace/zf_deseq_ts1/te_de/FAILED 2>/dev/null'
echo "== thinkstation2: analysis"
r thinkstation2 'echo "watcher alive: $(pgrep -f "[z]f_analyze_ts2" | wc -l)"; tail -1 ~/workspace/Geneformer_TE/runs/STATS_STATUS; ls ~/workspace/Geneformer_TE/runs/zf_validation_20261001/figures 2>/dev/null'
exit 0
