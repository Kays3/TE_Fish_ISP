# Where each part came from

Source: `Kays3/Geneformer_TE` (private), `origin/main` at `0f0d812085258199e6b2df11c907ff2aaba02054`, copied on
25 October 2026. Nothing was deleted or changed there. Every migrated file is byte-identical to that commit and keeps
its relative path. The last commit that touched each source path on `main` is given.

| Here | Source path | Last commit on main |
|---|---|---|
| `scripts/05_bulk_linear_isp.py` | same | `fb700ab` |
| `scripts/05b_te_axis_diagnostics.py`, `05c_compare_random_controls.py` | same | `a487ca9` |
| `scripts/06_te_axis_treatment.py` | same | `d8a91aa` |
| `scripts/07_treatment_gene_axis.py` | same | `71654f0` |
| `scripts/08_zebrafish_validation.py` | same | `c471d8c` |
| `scripts/08b_gse234993_sensitivity.py` | same | `17c21b8` |
| `scripts/08c` to `08f`, `04_bulk_validation.py`, `04_bulk_te_deseq.R`, `04b_bulk_te_deseq_mapcov.R` | same names | earlier |
| `scripts/utils/` | same | `6d8a752` |
| `config/` | same | `665fe99` |
| `pipelines/zebrafish_validation/` (16 files) | same | `73d46d0` |
| `docs/runs/2026-10-01_fish_bulk_isp/` | same | `2ad71cd` |
| `docs/runs/2026-10-01_te_axis_treatment/` | same | `2ad71cd` |
| `docs/runs/2026-10-01_gene_axis/` | same | `3c882d2` |
| `docs/runs/2026-10-01_zebrafish_validation/` | same | `74cc729` |
| `docs/BULK_LINEAR_ISP.md`, `docs/BULK_PERTURBATION_VALIDATION.md` | same | `526cc05` |
| `requirements.txt` | same | |
| `docs/context_from_source.md` | excerpts of `README.md`, `docs/DESIGN_CHANGES.md`, `docs/data-locations.md` | host names replaced |
| `data/small_inputs/` | four small TSV inputs of 05b and 06 as used on 1 October (sample sheet, two QC tables, the OA differential-expression table); not in git there | |
| `manifests/inputs.sha256.tsv` | verified transfer manifest of the archive that holds the large inputs | `manifests/README.md` |
| `docs/history.md` | the fish ISP history page of 24 October 2026 | |

Not new analysis: no table was recomputed, and the 05 run drew no figure, so none exists here.

Checked on 25 October 2026: all 114 migrated files (116 files compared by path; the other two, `README.md` and `.gitignore`, are this repository's own) are byte-identical to `origin/main` at `0f0d812`.
