#!/usr/bin/env bash
set -euo pipefail

# STAR's BAMsort temp files fan out to (threads x sort bins) open fds, which
# can exceed the default 1024 soft limit; raise it toward the hard limit
# before launching STAR.
ulimit -n 65536

if [[ $# -ne 5 ]]; then
  echo "Usage: $0 SAMPLES_TSV STAR_INDEX GENE_GTF TE_GTF OUTPUT_DIRECTORY" >&2
  exit 2
fi

manifest=$(realpath "$1")
star_index=$(realpath "$2")
gene_gtf=$(realpath "$3")
te_gtf=$(realpath "$4")
output_dir=$(realpath -m "$5")
threads=${THREADS:-14}
parallel_jobs=${PARALLEL_JOBS:-2}
tecount_jobs=${TECOUNT_JOBS:-3}
# Geneformer_TE zebrafish copy: an override file lets a running chain raise
# TEcount concurrency between waves without editing this script.
if [[ -s "$output_dir/TECOUNT_JOBS" ]]; then tecount_jobs=$(cat "$output_dir/TECOUNT_JOBS"); fi

mkdir -p "$output_dir"/{trimmed,bam,counts,qc}
pipeline_lock="$output_dir/.align-and-count.lock"
mkdir "$pipeline_lock" 2>/dev/null || {
  echo "Another align-and-count process already owns this output directory: $output_dir" >&2
  echo "If the previous process crashed, confirm it is stopped before removing $pipeline_lock" >&2
  exit 1
}
tecount_slots="$output_dir/counts/.tecount-slots"
rm -rf "$tecount_slots"
mkdir -p "$tecount_slots"
# Stale failure markers from a prior run would double-count against
# sample_count in the completion-wait loop below and could make it declare
# victory before every sample in *this* run has actually finished.
rm -f "$output_dir"/counts/*.tecount.failed

star_loaded=0
cleanup() {
  if (( star_loaded == 1 )); then
    STAR --genomeDir "$star_index" --genomeLoad Remove --runThreadN 1 >/dev/null 2>&1 || true
  fi
  rm -rf "$tecount_slots"
  rmdir "$pipeline_lock" 2>/dev/null || true
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

# Load the genome into shared memory once; concurrent per-sample STAR runs
# below reuse it via --genomeLoad LoadAndKeep instead of each paying the
# reload cost and duplicating it in RAM.
STAR --genomeDir "$star_index" --genomeLoad LoadAndExit --runThreadN "$threads"
star_loaded=1

process_sample() {
  local sample_id=$1 read1=$2 read2=$3 strandedness=$4
  local complete_marker="$output_dir/counts/${sample_id}.tecount.complete"
  # On a first run this per-sample directory does not exist yet. Guard the
  # find rather than letting it fail: under `set -e` with `pipefail`, a failed
  # find propagates through the pipe and aborts the whole run on sample one,
  # because `local` is declared separately so the assignment's exit status is
  # not masked.
  local existing_count_count=0
  if [[ -d "$output_dir/counts/$sample_id" ]]; then
    existing_count_count=$(find "$output_dir/counts/$sample_id" -maxdepth 1 \
      -type f -name '*.cntTable' -size +0c 2>/dev/null | wc -l)
  fi
  if [[ -s "$complete_marker" && "$existing_count_count" -eq 1 ]]; then
    echo "Skipping completed sample: $sample_id"
    return
  fi
  rm -f "$complete_marker"

  local fastp_args=(
    --in1 "$read1"
    --out1 "$output_dir/trimmed/${sample_id}_R1.fastq.gz"
    --thread "$threads"
    --html "$output_dir/qc/${sample_id}.fastp.html"
    --json "$output_dir/qc/${sample_id}.fastp.json"
  )
  local star_reads=("$output_dir/trimmed/${sample_id}_R1.fastq.gz")

  if [[ "$read2" != "-" ]]; then
    fastp_args+=(
      --in2 "$read2"
      --out2 "$output_dir/trimmed/${sample_id}_R2.fastq.gz"
      --detect_adapter_for_pe
    )
    star_reads+=("$output_dir/trimmed/${sample_id}_R2.fastq.gz")
  fi

  local bam="$output_dir/bam/${sample_id}.Aligned.sortedByCoord.out.bam"
  # Geneformer_TE zebrafish copy: reuse a finished alignment after a restart.
  if [[ -s "$bam.bai" && -s "$output_dir/bam/${sample_id}.Log.final.out" ]]; then
    echo "Reusing finished alignment: $sample_id"
  else
  fastp "${fastp_args[@]}"

  STAR \
    --genomeLoad LoadAndKeep \
    --limitBAMsortRAM "${BAM_SORT_RAM:-20000000000}" \
    --runThreadN "$threads" \
    --genomeDir "$star_index" \
    --readFilesIn "${star_reads[@]}" \
    --readFilesCommand zcat \
    --outFileNamePrefix "$output_dir/bam/${sample_id}." \
    --outFilterMultimapNmax 100 \
    --winAnchorMultimapNmax 100 \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI AS NM MD

  samtools index -@ "$threads" "$bam"
  fi

  mkdir -p "$output_dir/counts/$sample_id"
  # Run TEcount in the background (and disown it) so this slot can start the
  # next sample's fastp/STAR instead of blocking on it. TEcount itself is
  # capped to tecount_jobs concurrent runs via mkdir-as-lock slot directories,
  # since it's slower than expected and left uncapped it piles up faster than
  # it drains, oversubscribing the CPUs already committed to fastp/STAR above.
  (
    slot=""
    while [[ -z "$slot" ]]; do
      for s in $(seq 0 $((tecount_jobs - 1))); do
        if mkdir "$tecount_slots/$s" 2>/dev/null; then
          slot="$s"
          break
        fi
      done
      [[ -n "$slot" ]] || sleep 5
    done
    # Testing the exit status (rather than a bare call) keeps set -e, which
    # this subshell inherits from the parent script, from skipping the rmdir
    # below on failure — that would permanently leak the slot and eventually
    # deadlock every other sample.
    if ! TEcount \
      --sortByPos \
      --format BAM \
      --mode multi \
      --stranded "$strandedness" \
      --BAM "$bam" \
      --GTF "$gene_gtf" \
      --TE "$te_gtf" \
      --project "$sample_id" \
      --outdir "$output_dir/counts/$sample_id" \
      > "$output_dir/counts/${sample_id}.tecount.log" 2>&1; then
      touch "$output_dir/counts/${sample_id}.tecount.failed"
    else
      count_count=$(find "$output_dir/counts/$sample_id" -maxdepth 1 \
        -type f -name '*.cntTable' -size +0c | wc -l)
      if [[ "$count_count" -ne 1 ]]; then
        echo "Expected exactly one non-empty cntTable, found $count_count" \
          >> "$output_dir/counts/${sample_id}.tecount.log"
        touch "$output_dir/counts/${sample_id}.tecount.failed"
      else
        count_file=$(find "$output_dir/counts/$sample_id" -maxdepth 1 \
          -type f -name '*.cntTable' -size +0c)
        marker_tmp="$output_dir/counts/${sample_id}.tecount.complete.tmp.$$"
        {
          printf 'sample_id=%s\n' "$sample_id"
          printf 'completed_at=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
          printf 'count_file=%s\n' "$count_file"
          sha256sum "$count_file"
        } > "$marker_tmp"
        mv "$marker_tmp" "$complete_marker"
      fi
    fi
    rmdir "$tecount_slots/$slot"
  ) &
  disown
}

sample_count=$(tail -n +2 "$manifest" | cut -f1 | grep -c .)
running=0
tail -n +2 "$manifest" |
while IFS=$'\t' read -r sample_id read1 read2 treatment phenotype strandedness _; do
  [[ -n "$sample_id" ]] || continue
  read1=$(realpath "$read1")
  [[ "$read2" == "-" ]] || read2=$(realpath "$read2")

  process_sample "$sample_id" "$read1" "$read2" "$strandedness" &

  running=$((running + 1))
  if (( running >= parallel_jobs )); then
    wait -n
    running=$((running - 1))
  fi
done
wait  # drain remaining fastp/STAR/index slots (not the backgrounded TEcount jobs)

echo "Waiting for background TEcount jobs to finish writing counts..."
while :; do
  finished=$(find "$output_dir/counts" -maxdepth 1 \
    -name '*.tecount.complete' -type f -size +0c 2>/dev/null | wc -l)
  failed=$(find "$output_dir/counts" -maxdepth 1 -name '*.tecount.failed' 2>/dev/null | wc -l)
  (( finished + failed >= sample_count )) && break
  sleep 30
done
if (( failed > 0 )); then
  echo "WARNING: TEcount failed for $failed sample(s), continuing with the rest:" >&2
  find "$output_dir/counts" -maxdepth 1 -name '*.tecount.failed' 2>/dev/null | sed 's#.*/##;s/\.tecount\.failed$//' >&2
fi

multiqc "$output_dir" --outdir "$output_dir/qc/multiqc" --force
