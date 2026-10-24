# Run record: zebrafish ground truth for TE-regulator perturbation (option 2)

**Status at 2026-10-02 00:30 UTC:** complete. 112/112 counted (`COUNTS_READY`
21:00 UTC), DESeq2 17 studies 0 failed (`DESEQ_DONE` 21:07 UTC), 08 exit 0, rerun
twice after figure fixes. Read-1-only (GPU host, 19:13 UTC) and DC-18 mapping-covariate
runs done. `REPORT.md` written; awaiting review. Approved by the user: full panel, about 281 GB.
(Earlier status, 2026-10-01 10:10 UTC: downloads about 3% done.)

## Design

- **Panel.** 17 contrasts from 14 GEO series, 112 runs
  (`config/zebrafish_perturbation_studies.tsv`), found by a GEO query for
  zebrafish RNA-seq naming the configured regulators. Each study was then
  curated by hand: genetic loss or gain, bulk RNA-seq, replicated arms.
  - Covered: uhrf1 KO ×4, dnmt1 KO ×5 (liver, larva, and EC, HE and HSC
    sorted cells), dnmt3 morphants ×2, kdm1a/lsd1 KO ×2 alleles, ezh2
    maternal-zygotic KO, mettl3 morphant, atrx KO (MPNST), and
    human-UHRF1 overexpression in liver.
  - Excluded: scRNA-seq, drug-only, methylation-only and unreplicated
    studies (GSE103913, GSE93927).
- **References.** Ensembl GRCz11 release 112 primary assembly and GTF
  (MD5-checked against Ensembl CHECKSUMS). The TE GTF comes from UCSC
  danRer11 rmsk via `make_te_gtf.py`: 3,562,104 loci in 11,255 subfamilies.
  Alt-haplotype loci are dropped, and `chrUn_*vN` is mapped to Ensembl `*.N`.
- **Quantification.** The TE_Ocean_Acidification pipeline (fastp, STAR
  2.7.11b, TEcount 2.2.4 `--mode multi`), unstranded for every study.
- **DE.** DESeq2 1.46.0 with apeglm 1.28.0, from the frozen OA environment on
  compute host 1, design `~ role`.
- **Analysis.** `scripts/08_zebrafish_validation.py` reports, per study:
  - TE share of reads (Welch t, 95% CI);
  - TE families significantly up vs down (binomial test);
  - class-level medians (Wilcoxon);
  - the perturbed gene's own fold change.

  Across studies, it tests concordance with two predictions: "silencer loss
  derepresses TEs", and the fish bulk-ISP direction from 05c.

## Where things are

| Item | Location |
|---|---|
| Downloads, reference, BAMs, counts | alignment host `$ZF_WORK/` (`fastq/`, `reference/`, `rnaseq/`, `logs/STATUS`) |
| Pipeline env (x86 copy of the OA env) | alignment host `$ZF_WORK/env/te-ocean-acidification` |
| Shared outputs | NAS `${DATA_ROOT}/te-oa-full/zf-te-validation-20261001/` (`counts/`, `qc/`, `te_de/`, `COUNTS_READY`, `DESEQ_DONE`) |
| DESeq2 working dir | compute host 1 `$ZF_DESEQ/` |
| Analysis outputs | compute host 2 `$REPO_ROOT/runs/zf_validation_20261001/` |

## Watchers (each waits for the previous marker)

| Host | Process | Waits for | Writes |
|---|---|---|---|
| alignment host | `zf_download.sh` | — | `logs/STATUS: download exit=…` |
| alignment host | `zf_align_waves.sh` (replaced `zf_align.sh` at 12:31 UTC) | `reference done`; aligns every verified run in 10-min waves while downloads continue | `logs/align_waves.log`; `logs/STATUS: align exit=…` once downloads end and all verified runs are counted |
| alignment host | `zf_postalign.sh` | `align exit` | NAS `COUNTS_READY` |
| compute host 1 | `zf_deseq_ts1.sh` | NAS `COUNTS_READY` | NAS `DESEQ_DONE` |
| compute host 2 | `zf_analyze_ts2.sh` | NAS `DESEQ_DONE` | `runs/STATS_STATUS: 08 exit=…` |

## Resume after a disconnect

17:26 to 17:35 UTC, changes after the 17:30 pickup:

- **DESeq2 thread cap (operational, not a design change).** `zf_deseq_ts1.sh`
  caps BLAS/OpenMP at 4 threads (`DESEQ_THREADS`, commit `ce5d796`) on the
  shared host; `scripts/04_bulk_te_deseq.R` uses no BiocParallel. The
  compute host 1 watcher was restarted at 17:25 UTC: the old process group was
  stopped by PGID, and the new watcher (pid 2816374) started while
  `COUNTS_READY` was still absent.
- **Read-1-only realignment on GPU host.** Moved from alignment host to GPU host (aarch64, 12
  threads, NAS mounted) in a new aarch64 env with the alignment host tool
  versions; STAR index copied from alignment host and sha256-checked. Outputs go to NAS
  `zf-te-validation-20261001/r1only/{counts,qc}`, then `r1only/R1ONLY_READY`
  (`r1only exit=<rc> nas=<n>`, 9 libraries expected). Env path, installed
  versions, start/finish and CPU-h to be added when Oscar reports. After
  `R1ONLY_READY`, DESeq2 on the r1only counts runs by hand on compute host 1
  with the same clutch-blocked sheet (`~ block + role`).
- **Control SRR24921923** (GSE234993, clutch 1): not yet aligned at 17:35 UTC; the mean and
  tests above are to be updated once it is.

17:30 UTC (2026-10-01, host clock), pickup by a new session:

- Chain healthy: 220/220 FASTQ verified, 52/112 counted, 0 failed, all
  watchers alive. The current wave's manifest holds 74 runs (written at 14:07
  when 74 were verified); the next wave takes all 112.
- **GSE234993 mapping check (HANDOFF §4).** STAR unique-mapping %, 8 of 9
  libraries aligned (control SRR24921923, clutch1, still pending):

  | Run | Role | Clutch | Unique % | Multi % | Too short % |
  |---|---|---|---|---|---|
  | SRR24921922 | control | 2 | 22.62 | 20.83 | 56.18 |
  | SRR24921921 | control | 3 | 30.07 | 28.81 | 40.82 |
  | SRR24921920 | control | 4 | 31.71 | 25.27 | 42.71 |
  | SRR24921918 | perturbed | 2 | 23.13 | 24.93 | 51.56 |
  | SRR24921917 | perturbed | 3 | 21.87 | 16.30 | 61.39 |
  | SRR24921909 | perturbed | 4 | 19.26 | 17.02 | 63.29 |
  | SRR24921910 | perturbed | 5 | 19.36 | 20.09 | 60.15 |
  | SRR24921919 | perturbed | 1 | 20.99 | 21.65 | 56.95 |

  Mutant mean 20.9% vs control 28.1% (difference -7.2 points; Welch
  p = 0.11, Mann-Whitney p = 0.071, n = 5 vs 3). Within clutches 3 and 4 the
  mutant maps 8 to 12 points lower; clutch 2 is equal. Not significant at this
  n, but large relative to the ~25% baseline and driven by the "too short"
  fraction, so the read-1-only sensitivity realignment is queued.
- **Superseded at 17:26 UTC: the read-1-only realignment moved to GPU host** (see the 17:26 entry above). The alignment host watcher was stopped before it started,
  and alignment host `logs/STATUS` reads 'r1only HELD: moved to GPU host'. `zf_r1only.sh`
  stays in the repo as the method recipe, but it is not the path that ran.
- **`zf_r1only.sh` (alignment host, detached; stopped 17:26 UTC).** Waits for `align exit` (the main
  `align-and-count.sh` removes STAR's shared genome on exit, so the two must
  not overlap), realigns the 9 GSE234993 runs single-end from read 1 into
  `rnaseq_r1only/`, copies counts and STAR logs to NAS `r1only/`, and writes
  `R1ONLY_READY` and `logs/STATUS: r1only exit=…`. DESeq2 on these counts is
  run by hand afterwards with the same clutch-blocked sheet.

16:15 UTC, interim downstream test on the two fully counted studies (scratch
area; the NAS markers were not touched):

- **dnmt1 KO liver (3 vs 3).** 112 TE families up vs 26 down (binomial
  p < 1e-15), led by LTR families (median log2FC +0.14; 89% of significant
  LTR families up). TE share +16% [-15, +56], p = 0.22, underpowered.
- **uhrf1 KO liver (5 vs 5) under `~ role`.** 3 DE genes and 0 TE families.
  The study pairs each mutant with a sibling from the same clutch, and the
  mutant TE share is higher in 4 of 5 clutches. Under `~ block + role`:
  1,484 DE genes and 29 TE families.
- **Fix.** `config/zebrafish_perturbation_studies.tsv` has a `block` column
  (clutch) for the three clutch-paired studies (uhrf1 liver 2020 and 2023,
  hUHRF1 OE liver). `scripts/04_bulk_te_deseq.R` uses `~ block + role` when a
  block spans both roles. The compute host 1 watcher was restarted with the
  block-aware sheet builder.
- **Perturbed-gene check.** Neither uhrf1 (insertional hi272) nor dnmt1
  (point mutant s904) mRNA falls (log2FC +0.08 and -0.31, both NS). This is
  expected for these alleles, so the check is uninformative for them.


14:10 UTC restart:

- The first TEcount runs took 30 to 60 min each, and with 6 slots TEcount
  was the bottleneck (about 14 h for 112 samples).
- `pipelines/zebrafish_validation/align-and-count.sh` is a copy of the TE_OA
  script with two additions, shown by `diff` against the original:
  1. TEcount concurrency is read from `rnaseq/TECOUNT_JOBS`, now set to 16.
  2. A finished alignment (BAM index plus `Log.final.out`) is reused instead
     of being recomputed.
- The wave loop was restarted. 5 counted samples and 25 finished alignments
  were kept; 2 half-aligned samples were reset.
- Mapping (STAR unique) was 83 to 92% for most libraries; three were 62 to
  66%, to check in QC. The first count table had TE = 7.0% of counted reads.


Pickup at 2026-10-01 12:31 UTC: 81 of 220 files verified (100 GB,
about 10 MB/s, no failures). Alignment was switched to waves so it overlaps
with the remaining downloads (about 5 h left). Wave 1 started on 40 runs.


**Status** (the polling script `zf_status.sh` of the source repository is not copied here). It reported downloads
verified out of 220, samples counted out of 112, the NAS markers, the
DESeq2 tables and the analysis status, and whether each watcher is still
alive.

Snapshot at 2026-10-01 10:10 UTC, when the user disconnected:

- 13 GB downloaded, 0 of 220 files verified. Files are 1–4 GB at about
  0.6 MB/s per stream, so the first checksums land after roughly an hour.
- Reference done; 2.0 TB free on alignment host.
- All watchers alive.


1. Check progress:
   - `ssh alignment host 'cat $ZF_WORK/logs/STATUS; ls $ZF_WORK/fastq/*.ok | wc -l'` (220 files expected)
   - `ls ${DATA_ROOT}/te-oa-full/zf-te-validation-20261001/`
2. If a watcher died, rerun its script. Each one is idempotent:
   - Downloads resume, and verified files carry a `.ok` marker.
   - `align-and-count.sh` skips samples with a `.tecount.complete` marker.
3. When `08 exit=0` appears on compute host 2, pull
   `runs/zf_validation_20261001/{tables,figures}` into this folder and
   write the report.

## Known caveats to carry into the report

- Morpholino knockdowns (dnmt3, mettl3) that block translation need not
  lower target mRNA.
- The human UHRF1 transgene is not a zebrafish gene, so the
  perturbed-gene check does not apply to the OE study.
- The atrx study compares MPNST tumours on a tp53-mutant background, so the
  contrast includes tumour heterogeneity.
- Developmental stages and tissues differ across studies. Results are
  per study, never pooled into one model.

## Data moved

- **2026-10-04 12:27 UTC.** The whole working folder `$ZF_WORK/` on the lab compute server (section "Where things are": downloads, reference, BAMs, counts, pipeline env) was copied to `$SERVER_ARCHIVE/workspace/zf_te_validation/` on workstation D and verified there: 1,992 files, 834,631,865,258 bytes. The shared outputs under `te-oa-full/` on the lab NAS did not move. The old copy stays on the lab compute server until this record is merged. `docs/data-locations.md` holds the full record.
- **2026-10-05 06:08 UTC.** The old copy of `$ZF_WORK/` on the lab compute server was removed on the author's word (1,988 files, 830,040,140,665 bytes; the pipeline environment, 92,238 files, is held as `env.tar` and `.venv.tar`). The `ssh alignment host` status command under "Check progress" above no longer works. The folder now lives at `$SERVER_ARCHIVE/workspace/zf_te_validation/` on workstation D. The pipeline environment is x86_64 and does not run on workstation D (aarch64). `docs/data-locations.md` holds the record.
