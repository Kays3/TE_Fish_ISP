# Zebrafish ground-truth pipeline (option 2)

17 contrasts from 12 GEO series (112 runs) in which a TE regulator is knocked
out, knocked down or overexpressed in zebrafish. They are used as ground
truth for the fish bulk ISP. The study list with run accessions, roles and
sizes is `config/zebrafish_perturbation_studies.tsv`; `zf_download_list.tsv`
lists every FASTQ with its ENA MD5.

| Step | Host | Script |
|---|---|---|
| Download FASTQs (MD5-verified, resumable) and references | pSS | `zf_download.sh` |
| GRCz11 STAR index; TE GTF from UCSC danRer11 RepeatMasker (`make_te_gtf.py`; UCSC to Ensembl contig names, alt haplotypes dropped) | pSS | `zf_build_reference.sh` |
| fastp, STAR and TEcount `--mode multi`, unstranded (TE_Ocean_Acidification `align-and-count.sh`) | pSS | `zf_align.sh` |
| Publish count tables, alignment QC and gene names to the NAS | pSS | `zf_postalign.sh` |
| Per-study DESeq2 with apeglm (frozen OA environment: DESeq2 1.46.0, apeglm 1.28.0) | thinkstation1 | `zf_deseq_ts1.sh` |
| Statistics and figures (`scripts/08_zebrafish_validation.py`) | thinkstation2 | `zf_analyze_ts2.sh` |

Each step waits for the previous step's marker, so the chain runs
unattended. Shared outputs go to `/mnt/nas/kaisar/te-oa-full/zf-te-validation-20261001/`.

Counting is unstranded for every study, because library strandedness is not
reliably annotated across these submissions. This costs some resolution for
TEs antisense to genes.
