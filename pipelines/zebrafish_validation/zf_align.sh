#!/usr/bin/env bash
# Wait for downloads + reference, then align and count all zebrafish runs with the
# TE_Ocean_Acidification pipeline (fastp -> STAR -> TEcount --mode multi).
set -uo pipefail
W=$ZF_WORK; E=$W/env/te-ocean-acidification; R=$W/reference
mkdir -p $W/wrappers
cat > $W/wrappers/TEcount <<WR
#!/usr/bin/env bash
exec $E/bin/python $E/bin/TEcount "\$@"
WR
chmod +x $W/wrappers/TEcount
if [ ! -x $E/bin/multiqc ]; then printf '#!/usr/bin/env bash\necho "multiqc not installed; skipped"\n' > $W/wrappers/multiqc; chmod +x $W/wrappers/multiqc; fi
export PATH=$W/wrappers:$E/bin:$PATH
until grep -q "reference done" $W/logs/STATUS 2>/dev/null && grep -q "download exit" $W/logs/STATUS 2>/dev/null; do sleep 120; done
# sample sheet: only runs whose every FASTQ verified
$E/bin/python - <<PY
import csv, os
W = os.path.expanduser("$ZF_WORK")
rows = list(csv.DictReader(open(f"{W}/zebrafish_perturbation_studies.tsv"), delimiter="\t"))
out, missing = [], []
for r in rows:
    files = [f"{W}/fastq/" + os.path.basename(u) for u in r["ftp"].split(";")]
    if not all(os.path.exists(f + ".ok") for f in files):
        missing.append(r["run_accession"]); continue
    r1 = [f for f in files if f.endswith("_1.fastq.gz")] or files[:1]
    r2 = [f for f in files if f.endswith("_2.fastq.gz")]
    paired = r["layout"].upper() == "PAIRED" and r2
    out.append([r["run_accession"], r1[0], r2[0] if paired else "-", r["role"], r["study_id"], "no"])
with open(f"{W}/samples_all.tsv", "w") as f:
    f.write("sample_id\tread1\tread2\ttreatment\tparental_phenotype\tstrandedness\n")
    for o in out: f.write("\t".join(o) + "\n")
print(f"{len(out)} runs ready; missing: {missing}")
PY
THREADS=8 PARALLEL_JOBS=3 TECOUNT_JOBS=6 BAM_SORT_RAM=30000000000 \
  bash $W/align-and-count.sh $W/samples_all.tsv $R/star_index $R/genes.gtf $R/GRCz11_rmsk_TE.gtf $W/rnaseq > $W/logs/align.log 2>&1
echo "align exit=$? complete=$(ls $W/rnaseq/counts/*.tecount.complete 2>/dev/null | wc -l)" >> $W/logs/STATUS
