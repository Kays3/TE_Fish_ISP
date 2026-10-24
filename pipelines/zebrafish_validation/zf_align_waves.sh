#!/usr/bin/env bash
# Align and count in waves while downloads continue: each wave processes every
# run whose FASTQs are MD5-verified (align-and-count.sh skips completed samples).
# Finishes with "align exit=..." in logs/STATUS once downloads are done and every
# verified run has been processed.
set -uo pipefail
W=~/workspace/zf_te_validation; E=$W/env/te-ocean-acidification; R=$W/reference
mkdir -p $W/wrappers
cat > $W/wrappers/TEcount <<WR
#!/usr/bin/env bash
exec $E/bin/python $E/bin/TEcount "\$@"
WR
chmod +x $W/wrappers/TEcount
export PATH=$W/wrappers:$E/bin:$PATH
until grep -q "reference done" $W/logs/STATUS 2>/dev/null; do sleep 120; done
wave=0
while :; do
  dl_done=0; grep -q "download exit" $W/logs/STATUS 2>/dev/null && dl_done=1
  $E/bin/python - <<PY
import csv, os
W = os.path.expanduser("~/workspace/zf_te_validation")
rows = list(csv.DictReader(open(f"{W}/zebrafish_perturbation_studies.tsv"), delimiter="\t"))
out = []
for r in rows:
    files = [f"{W}/fastq/" + os.path.basename(u) for u in r["ftp"].split(";")]
    if not all(os.path.exists(f + ".ok") for f in files):
        continue
    r1 = [f for f in files if f.endswith("_1.fastq.gz")] or files[:1]
    r2 = [f for f in files if f.endswith("_2.fastq.gz")]
    paired = r["layout"].upper() == "PAIRED" and r2
    out.append([r["run_accession"], r1[0], r2[0] if paired else "-", r["role"], r["study_id"], "no"])
with open(f"{W}/samples_wave.tsv", "w") as f:
    f.write("sample_id\tread1\tread2\ttreatment\tparental_phenotype\tstrandedness\n")
    for o in out: f.write("\t".join(o) + "\n")
PY
  n_ready=$(($(wc -l < $W/samples_wave.tsv) - 1))
  n_done=$(ls $W/rnaseq/counts/*.tecount.complete 2>/dev/null | wc -l)
  if (( n_ready > n_done )); then
    wave=$((wave + 1))
    echo "$(date -u +%FT%TZ) wave $wave: ready=$n_ready done=$n_done" >> $W/logs/align_waves.log
    THREADS=8 PARALLEL_JOBS=3 TECOUNT_JOBS=6 BAM_SORT_RAM=30000000000 \
      bash $W/align-and-count.sh $W/samples_wave.tsv $R/star_index $R/genes.gtf $R/GRCz11_rmsk_TE.gtf $W/rnaseq >> $W/logs/align.log 2>&1
    echo "$(date -u +%FT%TZ) wave $wave exit=$?" >> $W/logs/align_waves.log
  elif (( dl_done == 1 )); then
    cp $W/samples_wave.tsv $W/samples_all.tsv
    echo "align exit=0 complete=$n_done ready=$n_ready" >> $W/logs/STATUS
    break
  else
    sleep 600
  fi
done
