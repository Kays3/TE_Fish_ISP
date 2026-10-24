# TE_Fish_ISP

In silico perturbation (ISP) of transposable-element (TE) regulators in fish, and the CellOracle-style network
propagation behind it, with the zebrafish ground-truth study that tested it. **Private repository: no public release,
no DOI.** Migrated on 25 October 2026 from `Kays3/Geneformer_TE`, where it was done on 1 to 13 October 2026; nothing
was deleted there.

## Question

Can natural variation in TE-silencer mRNA across unperturbed fish, propagated through a linear regulatory network
(a bulk adaptation of CellOracle, `scripts/05_bulk_linear_isp.py`), nominate the regulators of TE expression? The
fish bulk ISP was scored against 90 expression-matched random genes (ISP-STD-1 B1) and its predicted directions
against 17 zebrafish knockout, knockdown and overexpression contrasts.

## What was found (full account in `docs/history.md`)

| Study | Reading under ISP-STD-1 v1.1 |
|---|---|
| Fish bulk ISP, 30 regulators, 162 *A. polyacanthus* libraries | B1 not met: real regulators were not distinguishable from random genes, plain (393 vs 371 hits, p = 0.70) or with the global TE axis removed (62.5 vs 68, p = 0.17); 1 of 30 reached empirical p <= 0.05 (1.5 expected); 0 passed BH |
| Global TE-load axis | one axis carries 29.3% of within-stratum TE variance and tracks TE read fraction (r = 0.87); it drives the plain-mode directions |
| CO2 and the TE axis | developmental CO2 lowered the TE read share in gill by 11.4% [5.1, 17.2] (FDR 0.005); brain and liver inconclusive |
| Zebrafish ground truth, 17 contrasts, 112 runs | 0 of 16 loss-of-function studies shifted global TE share at FDR < 0.05; weakly powered, positive controls did not respond |
| Fish prediction against zebrafish response | not supported: Spearman rho = -0.18 (p = 0.50, n = 17) |

The data cannot show that these regulators do not act on fish TEs: unperturbed mRNA variation is a poor proxy for
silencing activity, and neither study had a responding positive control (`docs/history.md`, section 4).

## Data sources

| Data | Where it comes from | In this repository |
|---|---|---|
| 162 TEcount tables, *A. polyacanthus* OA cohort (unpublished) | `te-oa-full/counts-singlebatch/all/` | path and sha256 in `manifests/inputs.sha256.tsv`; files stay outside git |
| Annotation `genes.gtf` (GCF_021347895.1) | `te-oa-full/reference/` | same |
| Zebrafish counts, DESeq2 tables, alignment QC (14 GEO series, 112 runs) | `te-oa-full/zf-te-validation-20261001/` | same; accessions in `config/zebrafish_perturbation_studies.tsv` |
| Sample sheet, two QC tables, OA differential-expression table | derived from the cohort | `data/small_inputs/` (under 1 MB) |
| Raw reads (834.6 GB) | GEO and ENA, MD5 in `pipelines/zebrafish_validation/zf_download_list.tsv` | not copied; one archived copy exists |

`DATA_ROOT` is the directory that contains `te-oa-full/`. On 5 October the folders moved from the NAS to an archive
host; the manifest records the hashes verified at that copy. Nothing over 50 MB is in git, and no file was copied
through a laptop.

## How to rerun

```bash
scripts/run_all.sh                        # figures: all 11 rebuilt from results/ in seconds, no data needed
scripts/run_all.sh verify                 # checksums (repository files; add DATA_ROOT=... for the large inputs)
DATA_ROOT=/path scripts/run_all.sh compute   # recompute tables from the counts: CPU only, heavy, see below
```

- **Environment:** `environment/README.md` pins Python 3.11.15 with matplotlib 3.11.2. Figures are byte-reproducible
  only with that lock.
- **Figures:** `scripts/make_figures.py` draws each figure from the saved tables with the original plotting code.
  `figures/INDEX.md` lists every figure, its source table and its checksum.
- **Recompute:** the `compute` stage lists the commands of the 2026-10-01 run records (05, 05b, 05c, 06, 07, 08).
  It was not re-executed when the repository was made; run it on a compute host, not a laptop (step 05 needs about
  20 cores for 5 to 20 minutes per run). The alignment and DESeq2 steps for zebrafish are in
  `pipelines/zebrafish_validation/` and are not rerun by the driver.
- **No GPU** is used anywhere in this repository.

## Layout

```
scripts/                 analysis (05 to 08f), figure code (make_figures.py, fig05/06_plots.py), driver, checksums
pipelines/zebrafish_validation/   download, STAR + TEcount, DESeq2 wrappers (paths through $ZF_WORK, $ZF_DESEQ)
config/                  project config, TE-regulator gene sets, the 17-contrast zebrafish manifest
configs and run records  results/2026-10-01_*/   configs/, tables/, REPORT.md or RUN_RECORD.md of each run
data/small_inputs/       four small TSV inputs
manifests/               sha256 of the large inputs, repository files and figures
environment/             Python and R package pins
figures/                 rebuilt figures + INDEX.md
docs/                    history.md, source_map.md, background notes
planned/                 acute TE response and warming: planned, not yet run
```

## Provenance

`docs/source_map.md` gives the source commit of every file. Path and host strings were replaced by placeholders
(five provenance JSONs differ from the source in those strings only). Result tables are the committed ones, unchanged.
Wording of results follows `ISP-STD-1` (status vocabulary and "not estimable", never "no effect").

## Next work (planned, not yet run)

`planned/acute_te_response/` and `planned/warming/` hold the layout and rules for the next two lines of work. They
contain no analysis.
