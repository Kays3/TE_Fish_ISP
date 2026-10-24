# Where each part came from

Source: `Kays3/Geneformer_TE` (private), `origin/main` at `0f0d812085258199e6b2df11c907ff2aaba02054`, copied on
25 October 2026. Nothing was deleted from the source. The last commit that touched each source path on `main` is
given so a reader can find the version that was copied.

| Here | Source path | Last commit on main |
|---|---|---|
| `scripts/05_bulk_linear_isp.py` | `scripts/05_bulk_linear_isp.py` | `fb700ab` |
| `scripts/05b_te_axis_diagnostics.py`, `05c_compare_random_controls.py` | same names in `scripts/` | `a487ca9` |
| `scripts/06_te_axis_treatment.py` | `scripts/06_te_axis_treatment.py` | `d8a91aa` |
| `scripts/07_treatment_gene_axis.py` | `scripts/07_treatment_gene_axis.py` | `71654f0` |
| `scripts/08_zebrafish_validation.py` | `scripts/08_zebrafish_validation.py` | `c471d8c` |
| `scripts/08b` to `08f`, `04_bulk_validation.py`, `04_bulk_te_deseq.R`, `04b_bulk_te_deseq_mapcov.R` | same names in `scripts/` | `17c21b8` (08b) and earlier |
| `scripts/utils/` | `scripts/utils/` | `6d8a752` |
| `config/` | `config/` | `665fe99` |
| `pipelines/zebrafish_validation/` | `pipelines/zebrafish_validation/` (without `zf_status.sh`, which only polls hosts) | `73d46d0` |
| `results/2026-10-01_fish_bulk_isp/` | `docs/runs/2026-10-01_fish_bulk_isp/` (without figures) | `2ad71cd` |
| `results/2026-10-01_te_axis_treatment/` | `docs/runs/2026-10-01_te_axis_treatment/` | `2ad71cd` |
| `results/2026-10-01_gene_axis/` | `docs/runs/2026-10-01_gene_axis/` | `3c882d2` |
| `results/2026-10-01_zebrafish_validation/` | `docs/runs/2026-10-01_zebrafish_validation/` | `74cc729` |
| `docs/background/` | `docs/BULK_LINEAR_ISP.md`, `docs/BULK_PERTURBATION_VALIDATION.md` | `526cc05` |
| `data/small_inputs/` | the four small TSV inputs of steps 05b and 06 (sample sheet, two QC tables, the OA differential-expression table), as used on 1 October | not in git there |
| `manifests/inputs.sha256.tsv` | archive manifest of the large inputs | see `manifests/README.md` |
| `docs/history.md` | the fish ISP history page of 24 October 2026 | |

## Changes made while copying

1. **Paths and host names.** Home paths, the NAS mount and compute host names were replaced by placeholders in
   text files: `$REPO_ROOT`, `$ZF_WORK`, `$ZF_DESEQ`, `${RSCRIPT:-Rscript}`, `${DATA_ROOT}`, "compute host 1",
   "compute host 2", "alignment host", "GPU host". This also touches the `config_path` and `counts` strings inside
   the `*_provenance.json` files, so those five files differ from the source in those strings only.
2. **Access notes removed.** The "Access" paragraph of the fish run record and two sentences of the zebrafish run record that described how to reach the compute hosts were replaced by neutral text.
3. **Figures.** The source figures were drawn inline by the analysis scripts from raw inputs. Here
   `scripts/make_figures.py` redraws them from the saved tables, using the original plotting code (07 and 08 by
   import, 06 by a verbatim copy in `fig06_plots.py`). On the host that made the originals, with the pinned
   environment, the ten PNGs are byte-identical to the source PNGs (see `figures/INDEX.md`).
4. **One new figure**, `05_fig1_real_vs_random`, plots the saved 05c values. The 2026-10-01 run drew none.
5. **Not copied:** Geneformer-model steps (01 to 03, 09 to 11), the human T-cell runs, and everything about the
   SCLC, NSCLC and E3 work. The `config/config.yaml` copied here is the shared project file; only its `bulk_isp`,
   `bulk_validation` and path sections are used by the fish scripts.
