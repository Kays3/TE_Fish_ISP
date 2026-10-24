#!/usr/bin/env bash
# compute host 1: per-study DESeq2 (apeglm) in the frozen OA environment (DESeq2 1.46.0, apeglm 1.28.0).
set -uo pipefail
# Shared host: cap BLAS/OpenMP threads so a GPU job's data loading is not starved.
export OPENBLAS_NUM_THREADS=${DESEQ_THREADS:-4} OMP_NUM_THREADS=${DESEQ_THREADS:-4} MKL_NUM_THREADS=${DESEQ_THREADS:-4}
N=${DATA_ROOT}/te-oa-full/zf-te-validation-20261001
G=$REPO_ROOT
R=${RSCRIPT:-Rscript}
L=$ZF_DESEQ; mkdir -p $L/sheets $L/te_de $L/counts
until [ -s $N/COUNTS_READY ]; do sleep 300; done
cp $N/counts/*.cntTable $L/counts/
python3 - <<PY
import csv, collections
ready = {r["sample_id"] for r in csv.DictReader(open("$N/samples_all.tsv"), delimiter="\t")}
man = list(csv.DictReader(open("$G/config/zebrafish_perturbation_studies.tsv"), delimiter="\t"))
by = collections.defaultdict(list)
for r in man:
    if r["run_accession"] in ready:
        by[r["study_id"]].append(r)
for study, rs in by.items():
    with open(f"$L/sheets/{study}.tsv", "w") as f:
        f.write("sample_id\ttreatment\tblock\n")
        for r in rs: f.write(f"{r['run_accession']}\t{r['role']}\t{r.get('block', '')}\n")
print(len(by), "studies")
PY
for s in $L/sheets/*.tsv; do
  st=$(basename $s .tsv)
  $R $G/scripts/04_bulk_te_deseq.R $L/counts $s $L/te_de/$st.tsv > $L/te_de/$st.log 2>&1 || echo "FAILED $st" >> $L/te_de/FAILED
done
mkdir -p $N/te_de && cp $L/te_de/*.tsv $L/te_de/*.txt $N/te_de/ 2>/dev/null
echo "studies=$(ls $L/te_de/*.tsv | wc -l) failed=$(cat $L/te_de/FAILED 2>/dev/null | wc -l)" > $N/DESEQ_DONE
