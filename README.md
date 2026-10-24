# TE_Fish_ISP

In silico perturbation (ISP) of transposable-element (TE) regulators in fish, and the CellOracle-style network
propagation behind it, with the zebrafish ground-truth study that tested it. **Private repository: no public release,
no DOI.** Copied on 25 October 2026 from `Kays3/Geneformer_TE`, where the work was done on 1 to 13 October 2026.
**Source commit: `0f0d812085258199e6b2df11c907ff2aaba02054` (`origin/main`).** Nothing was deleted or changed there.

## Question

Can natural variation in TE-silencer mRNA across unperturbed fish, propagated through a linear regulatory network
(a bulk adaptation of CellOracle, `scripts/05_bulk_linear_isp.py`), nominate the regulators of TE expression? The
fish bulk ISP was scored against 90 expression-matched random genes (ISP-STD-1 B1) and its predicted directions
against 17 zebrafish knockout, knockdown and overexpression contrasts.

## What was found (full account in `docs/history.md`)

| Study | Reading under ISP-STD-1 v1.1 |
|---|---|
| Fish bulk ISP, 30 regulators, 162 *A. polyacanthus* libraries | B1 not met (regulators not distinguishable from random genes, and only 3 matched controls per regulator against the 20 per gene B1 asks for, so per-gene status is `not_estimable`): plain (393 vs 371 hits, p = 0.70) or with the global TE axis removed (62.5 vs 68, p = 0.17); 1 of 30 reached empirical p <= 0.05 (1.5 expected); 0 passed BH |
| Global TE-load axis | one axis carries 29.3% of within-stratum TE variance and tracks TE read fraction (r = 0.87); it drives the plain-mode directions |
| CO2 and the TE axis | developmental CO2 lowered the TE read share in gill by 11.4% [5.1, 17.2] (FDR 0.0056); brain and liver inconclusive |
| Zebrafish ground truth, 17 contrasts, 112 runs | 0 of 16 loss-of-function studies shifted global TE share at FDR < 0.05; weakly powered, positive controls did not respond |
| Fish prediction against zebrafish response | not supported: Spearman rho = -0.18 (p = 0.50, n = 17) |

The data cannot show that these regulators do not act on fish TEs: unperturbed mRNA variation is a poor proxy for
silencing activity, and neither study had a responding positive control (`docs/history.md`, section 4).

## Layout

The migrated files keep the relative paths they have in `Geneformer_TE`, because the scripts load each other and
`utils.common` by file path and the run records link to each other.

```
scripts/            04 to 08f analysis scripts, utils/; new here: make_figures.py, fig06_plots.py, run_all.sh, checksum scripts
pipelines/zebrafish_validation/   download, STAR + TEcount, DESeq2 wrappers (16 files, byte-identical)
config/             config.yaml (verbatim), gene_sets_apoly.yaml, zebrafish_perturbation_studies.tsv
docs/runs/2026-10-01_{fish_bulk_isp,te_axis_treatment,gene_axis,zebrafish_validation}/   run records, configs, tables, figures
docs/BULK_LINEAR_ISP.md, BULK_PERTURBATION_VALIDATION.md     method notes
docs/history.md, context_from_source.md, source_map.md        new here: what happened, copied context, provenance
data/small_inputs/  four small TSV inputs of 05b and 06 (new here; they were not in git)
manifests/          sha256 of the 568 large inputs, of repository files and of the committed figures
environment/        Python lock and DESeq2 package list; requirements.txt is the source file the lock is based on
figures/INDEX.md    the ten figures, their tables, code and checksums
planned/            acute TE response and warming: planned, not yet run
```

## Data sources

| Data | Where it comes from | In this repository |
|---|---|---|
| 162 TEcount tables, *A. polyacanthus* ocean-acidification cohort (unpublished) | `te-oa-full/counts-singlebatch/all/` | path and sha256 in `manifests/inputs.sha256.tsv`; stays outside git |
| Annotation `genes.gtf` (GCF_021347895.1) | `te-oa-full/reference/` | same |
| Zebrafish counts, DESeq2 tables, alignment QC (14 GEO series, 112 runs) | `te-oa-full/zf-te-validation-20261001/` | same; accessions in `config/zebrafish_perturbation_studies.tsv` |
| Sample sheet, two QC tables, OA differential-expression table | derived from the cohort | `data/small_inputs/` (under 1 MB) |
| Raw reads (834.6 GB) | GEO and ENA, MD5 in `pipelines/zebrafish_validation/zf_download_list.tsv` | not copied; one archived copy exists |

`DATA_ROOT` is the directory that contains `te-oa-full/`. The source repository's `docs/data-locations.md` says
`te-oa-full/` has not moved; it has. On 5 October 2026 the folders were copied to an archive host with a sha256
check and the NAS copies were removed (`docs/context_from_source.md`, last section). The manifest takes its hashes from
that verified copy. No file was read or copied through a laptop. The vendored `align-and-count.sh` is the
`TE_Ocean_Acidification` count pipeline as used on 2026-10-01; that repository's current file differs.

## How to rerun

```bash
scripts/run_all.sh                           # rebuild the 10 figures from docs/runs/*/tables and compare with the committed PNGs
scripts/run_all.sh verify                    # checksums (repository files; add DATA_ROOT=... for the large inputs)
DATA_ROOT=/path scripts/run_all.sh compute   # recompute tables from the counts: CPU only, heavy, see below
```

- **Environment:** `environment/README.md` pins Python 3.11.15 with matplotlib 3.11.2. Figures are byte-reproducible
  only with that lock.
- **Compute stage:** lists the commands of the 2026-10-01 run records (05, 05b, 05c, 06, 07, 08, 08c, 08d, 08f). It
  was not re-executed when this repository was made. Step 05 needs about 20 cores for 5 to 20 minutes per run. The
  zebrafish alignment and DESeq2 steps (`pipelines/zebrafish_validation/`) are not rerun by the driver; 08b and 08e
  need DESeq2 working folders that are not in the manifest.
- **No GPU** is used anywhere in this repository.

## Provenance notes

- **Byte identity.** Every migrated file is byte-identical to `Geneformer_TE` at `0f0d812`: the analysis scripts,
  `config/config.yaml`, the 16 pipeline files and every file under `docs/runs/2026-10-01_*` (including the committed
  figures). `config/config.yaml` keeps its single-cell `geneformer`, `classifier`, `isp` and `probe` sections; they are
  inert in this repository, and the provenance JSONs pin its sha256.
- **Home paths and host names.** Historical run records, provenance JSONs, run configs and pipeline scripts contain
  the absolute paths and host names of the original Linux run hosts. They are kept as they were run (they are data,
  not configuration of this repository). Every file this repository adds (this README, the driver, manifests, figure
  code, lock files, history, the data notes) has none; the one exception is `scripts/run_all.sh`, which names the historical NAS prefix as the `sed` pattern that rewrites it to `DATA_ROOT` when it renders configs.
- **Added files:** `scripts/make_figures.py`, `scripts/fig06_plots.py` (a verbatim copy of the inline plotting block
  of `06_te_axis_treatment.py`, which could not be called without the counts), `scripts/run_all.sh`,
  `scripts/verify_inputs.sh`, `scripts/write_checksums.sh`, `scripts/_sha.sh`, `data/small_inputs/`, `manifests/`,
  `environment/`, `docs/history.md`, `docs/context_from_source.md`, `docs/source_map.md`, `figures/INDEX.md`, `planned/`.

## Left behind in Geneformer_TE, and why

| What | Reason |
|---|---|
| `scripts/01` to `03`, `patches/geneformer-bf16-export.patch` | the single-cell Geneformer track: a template never run on real data, not specific to fish |
| `scripts/09` to `13`, `pipelines/lund_hnpc_te/`, `docs/runs/2026-10-02_*` | human datasets |
| `tests/` | all five test files cover scripts 09 to 13; **nothing in the repository tests the fish scripts**, and none were written here (out of scope for a migration) |
| `docs/graphical-abstract/` | covers the human runs too |
| `docs/HANDOFF.md`, `AGENT_PROMPT.md`, `DESIGN_CHANGES.md`, `data-locations.md` | repository-wide; the fish content is copied into `docs/context_from_source.md` |
| `notebooks/`, `results/` | empty placeholders |

## Next work (planned, not yet run)

`planned/acute_te_response/` and `planned/warming/` hold the layout and rules for the next two lines of work. They
contain no analysis.
