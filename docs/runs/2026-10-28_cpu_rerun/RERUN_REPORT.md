# CPU rerun of `scripts/run_all.sh compute` (28 October 2026)

## Summary

We reran the CPU analyses 05 to 08 from the archived counts, as the repository's `compute` stage describes, and compared every committed table with its rerun file byte for byte. **5 of the 35 recomputed tables are byte-identical** (a sixth file, `08e_clutch_deviation.tsv`, was copied in, not recomputed). Of the other 30, 22 differ only at floating-point rounding: the largest relative difference in any numeric cell is 9.2e-12. They print 16 to 17 significant digits, so a different summation order changes the last digits. 6 differ for known reasons: 5 are run metadata (paths, platform, commit, time), and `08_family_lfc.tsv.gz` has identical content and a different gzip header. The last two tables need a note:

- `06_phenotype_tests.tsv`: 18 of 108 rows change method and p: all six contrasts of three tissue-measure cells (gill `te_frac_strict`, gill `te_frac_all`, liver `axis_score`). 15 of the 18 also change the estimate by more than 1e-9; the three interaction rows keep it (difference about 1e-15). In these three cells the committed REML fit did not converge. Its family variance (5.0e-4, 1.0e-4 and 7.1e-6) is above the 1e-6 boundary threshold. The script's rule (variance below 1e-6, or no convergence) therefore replaced it with fish-level OLS (`ols_boundary`, a label that covers both cases). The rerun converged with the same family variance and kept REML (`mixedlm_reml`). The two methods use different reference distributions: `ols_boundary` takes p from t with the residual df, and `mixedlm_reml` from the normal (`scripts/06_te_axis_treatment.py` lines 140 to 145). Part of each p change therefore comes from the method, not only from the estimate. The smallest p change in the 18 rows is 0.0013. The largest estimate change is -0.0748 (committed) against -0.0710 (rerun), with p 0.162 against 0.197. Because the FDR family includes these rows, q also changes in 14 other rows (by up to 0.053) with no change of estimate or p; this gives the largest absolute difference in the table, 0.0905, in the `fdr` column of a switched row. In liver `te_frac_all` the family variance differs (4.3e-11 against 1.7e-10, near zero in both runs) in all six rows; two of them also differ in q. At a relative threshold of 1e-4, 36 of 108 rows differ: the 18 switched rows and 18 others. One row crosses p < 0.05 (gill `te_frac_strict`, `mean_of_treatments`, p 0.029 against 0.053), and no q < 0.05 call changes. The other two `ols_boundary` cells keep their method: brain `te_frac_strict` (no convergence in both runs) and liver `te_frac_all` (the only true variance-boundary case). The fit is fragile, and which method these rows use depends on the optimizer path. The rerun table is `06_phenotype_tests_rerun.tsv` in this folder. We did not check each claim in the run records that cites these rows.
- `08d_hUHRF1_OE_te_share_contributors.tsv`: the same values; three rows with the value 0.0 are in a different order.

`compute` did not run as committed. It failed twice, and we fixed three defects in `run_all.sh` in this change (section "Defects in run_all.sh"). We did not run the fixed script end to end. The rerun is the committed script up to step 05b, then the same commands with the fixes applied (section "What ran").

5 committed tables are not made by `compute` and were not checked (last rows of the table below). Steps 08b and 08e need the DESeq2 working folders and were not rerun.

The figures stage (`run_all.sh figures`, on the same machine and environment) gave 0 of 10 PNGs different from the committed figures. That stage draws from the committed tables, not from the rerun tables, so it tests the drawing code and the environment, not the computation.

## Setup

| Item | Value |
|---|---|
| Code | this repository at `1d882fb` (main), clean checkout |
| Machine | Linux x86_64 CPU server, 32 cores, `Linux-7.0.0-38-generic-x86_64-with-glibc2.39`; `THREADS=10`; no GPU |
| Python | 3.11.15, `environment/python-bulk.lock.txt` installed exactly (`uv pip freeze` equals the lock) |
| Inputs | the 568 files of `manifests/inputs.sha256.tsv` (538,394,492 bytes). They are no longer at their recorded NAS paths: they moved to an archive on 5 October 2026 (`MOVED_TO.txt` in the old folder). We copied them from the archive into the staging folder NAS `/mnt/nas/kaisar/fish-isp-rerun-20261028` (568 files, 538,394,492 bytes) and set `DATA_ROOT` to it. 568 of 568 match `manifests/inputs.sha256.tsv` (sha256 of that manifest file: `0fbf9a0d409683def5d2a8422b7c9a7855eaf4afeefc5224153f303b16c07537`). |
| Run times | 05 and 05b: 1 h 18 min 51 s, peak RSS 3.2 GB. 05c to 08d: 1 min 55 s, peak RSS 1.0 GB. 05b again with the QC fix: 30 s. 08f: 1 s. |

## How to repeat it

```bash
uv venv --python 3.11 .venv-bulk
uv pip install --python .venv-bulk/bin/python -r environment/python-bulk.lock.txt
DATA_ROOT=/mnt/nas/kaisar/fish-isp-rerun-20261028 scripts/verify_inputs.sh      # large inputs: all 568 OK
PYTHON=.venv-bulk/bin/python DATA_ROOT=/mnt/nas/kaisar/fish-isp-rerun-20261028 THREADS=10 scripts/run_all.sh compute
```

The staging folder stays until this change is merged; its removal after that is not yet decided. The lasting source of the 568 inputs is the archive copy named in `MOVED_TO.txt` under the old NAS folder `te-oa-full/counts-singlebatch/`; with that archive as `DATA_ROOT` the commands are the same.

## What ran

1. `run_all.sh compute`, 13:32 to 14:51 UTC. Steps 05 (all nine runs) and 05b completed; 05c failed on its input paths (defect 2).
2. The remaining commands of `stage_compute`, cut from `run_all.sh` at `1d882fb` with only the two 05c input lists changed to the three file pairs. 05c to 08d completed; 08f failed because `08e_clutch_deviation.tsv` was missing (defect 3).
3. 08f alone, with the committed `08e_clutch_deviation.tsv` copied into the rerun table folder first.
4. The comparison showed that `05b_axis_correlates.tsv` had other rows than the committed table (defect 1). 05b again with `--qc runs/qc_align.tsv`; its tables then matched to rounding. The comparison below uses this 05b output.

## Defects in run_all.sh (fixed here)

1. **05b used the wrong QC table.** `--qc runs/qc.tsv` gives the columns `raw_read_pairs`, `source_mapping_rate_pct` and `source_genes_expressed`. The committed `05b_axis_correlates.tsv` has `log_star_input_reads`, `unique_mapping_pct`, `multimapping_pct` and `unmapped_fraction_pct`, which are the columns of `qc_align.tsv`. Fixed to `--qc runs/qc_align.tsv`.
2. **05c inputs were a cartesian product.** `runs/fish_random{1,2,3}/tables/05_random{1,2,3}_bulk_isp_effects.tsv` expands in bash to 9 paths, 6 of which do not exist. The same holds for the `rmTEpc1` list. Fixed by writing the three pairs out.
3. **08f needs 08e, which `compute` does not make.** `08f_exploratory_concordance.py` reads `08e_clutch_deviation.tsv`. 08e needs the DESeq2 working folders, which this stage does not rerun. Fixed by copying the committed 08e table into the rerun folder when it is missing. 08f is then a recomputation from the committed 08e table, not from new 08e output.

## Result per output

The full machine-readable comparison is `comparison.tsv` (one row per committed file and rerun file of the same name). Paths in the first column are under `docs/runs/`. Rows where a committed file name matched a file of another run (the `no_outliers` copies of the 05 tables, and the copies of the 05c tables in `zf_validation_20261001/fish05c/`) are left out here; they are in `comparison.tsv`. Its last row, `figures/INDEX.md`, is not a table and is left out too.

| Committed output | Rerun output | Bytes | Numeric comparison | Note |
|---|---|---|---|---|
| `2026-10-01_fish_bulk_isp/tables/05_bulk_isp_provenance.json` | `runs/fish_bulk_isp_20261001/tables/05_bulk_isp_provenance.json` | differ | text cells differing: 7 | run metadata: config path and hash, counts path, platform, commit, time; spectral radius at rounding |
| `2026-10-01_fish_bulk_isp/tables/05_bulk_target_cv_r2.tsv` | `runs/fish_bulk_isp_20261001/tables/05_bulk_target_cv_r2.tsv` | differ | 1619 of 2235 cells differ; max abs 8.05e-16, max rel 5.83e-13 |  |
| `2026-10-01_fish_bulk_isp/tables/05_rmTEpc1_bulk_isp_provenance.json` | `runs/fish_bulk_isp_20261001_rmTEpc1/tables/05_rmTEpc1_bulk_isp_provenance.json` | differ | text cells differing: 7 | run metadata, as above |
| `2026-10-01_fish_bulk_isp/tables/05b_axis_correlates.tsv` | `runs/fish_bulk_isp_20261001/tables/05b_axis_correlates.tsv` | differ | 21 of 42 cells differ; max abs 9.99e-16, max rel 4.43e-13 |  |
| `2026-10-01_fish_bulk_isp/tables/05b_axis_summary.tsv` | `runs/fish_bulk_isp_20261001/tables/05b_axis_summary.tsv` | differ | 2 of 8 cells differ; max abs 3.33e-16, max rel 1.14e-15 |  |
| `2026-10-01_fish_bulk_isp/tables/05b_provenance.json` | `runs/fish_bulk_isp_20261001/tables/05b_provenance.json` | differ | text cells differing: 7 | run metadata, as above |
| `2026-10-01_fish_bulk_isp/tables/05b_regulator_axis_corr.tsv` | `runs/fish_bulk_isp_20261001/tables/05b_regulator_axis_corr.tsv` | differ | 94 of 210 cells differ; max abs 1.39e-15, max rel 4.06e-13 |  |
| `2026-10-01_fish_bulk_isp/tables/05c_plain_overall.tsv` | `runs/05c/05c_plain_overall.tsv` | differ | 1 of 28 cells differ; max abs 1.11e-16, max rel 8.35e-16 |  |
| `2026-10-01_fish_bulk_isp/tables/05c_plain_per_regulator.tsv` | `runs/05c/05c_plain_per_regulator.tsv` | differ | 46 of 600 cells differ; max abs 5e-16, max rel 2.5e-15 |  |
| `2026-10-01_fish_bulk_isp/tables/05c_rmTEpc1_overall.tsv` | `runs/05c/05c_rmTEpc1_overall.tsv` | differ | 1 of 28 cells differ; max abs 2.08e-16, max rel 2.48e-15 |  |
| `2026-10-01_fish_bulk_isp/tables/05c_rmTEpc1_per_regulator.tsv` | `runs/05c/05c_rmTEpc1_per_regulator.tsv` | differ | 58 of 600 cells differ; max abs 3.05e-16, max rel 5.07e-15 |  |
| `2026-10-01_fish_bulk_isp/tables/outlier_sensitivity_delete.tsv` | none | not produced | | not made by `compute` |
| `2026-10-01_gene_axis/tables/07_gene_axis_association.tsv` | `runs/fish_gene_axis_20261001/tables/07_gene_axis_association.tsv` | differ | 922 of 1736 cells differ; max abs 3.89e-15, max rel 4.17e-13 |  |
| `2026-10-01_gene_axis/tables/07_gene_mediation.tsv` | `runs/fish_gene_axis_20261001/tables/07_gene_mediation.tsv` | differ | 4752 of 6154 cells differ; max abs 4.53e-12, max rel 3.42e-13 |  |
| `2026-10-01_gene_axis/tables/07_gene_treatment_effects.tsv` | `runs/fish_gene_axis_20261001/tables/07_gene_treatment_effects.tsv` | differ | 2366 of 3038 cells differ; max abs 1.87e-14, max rel 8.39e-13 |  |
| `2026-10-01_gene_axis/tables/07_provenance.json` | `runs/fish_gene_axis_20261001/tables/07_provenance.json` | differ | text cells differing: 5 | run metadata, as above |
| `2026-10-01_gene_axis/tables/07_set_mediation.tsv` | `runs/fish_gene_axis_20261001/tables/07_set_mediation.tsv` | differ | 109 of 172 cells differ; max abs 1.74e-13, max rel 1.58e-12 |  |
| `2026-10-01_gene_axis/tables/07_set_treatment_shift.tsv` | `runs/fish_gene_axis_20261001/tables/07_set_treatment_shift.tsv` | differ | 127 of 312 cells differ; max abs 5.61e-15, max rel 4.54e-13 |  |
| `2026-10-01_te_axis_treatment/tables/06_de_alignment.tsv` | `runs/fish_te_axis_treatment_20261001/tables/06_de_alignment.tsv` | identical |  |  |
| `2026-10-01_te_axis_treatment/tables/06_de_alignment_with_shift.tsv` | `runs/fish_te_axis_treatment_20261001/tables/06_de_alignment_with_shift.tsv` | differ | 16 of 102 cells differ; max abs 1.89e-15, max rel 5.05e-15 |  |
| `2026-10-01_te_axis_treatment/tables/06_phenotype_tests.tsv` | `runs/fish_te_axis_treatment_20261001/tables/06_phenotype_tests.tsv` | differ | 633 of 918 cells differ; max abs 0.0905, max rel 3.03 | 18 of 108 rows change method and p (6 in each of gill te_frac_strict, gill te_frac_all and liver axis_score): the committed REML fit did not converge and the script used fish-level OLS; the rerun converged and kept REML. q changes in 14 more rows; the max abs diff is in `fdr` (see Summary) |
| `2026-10-01_te_axis_treatment/tables/06_provenance.json` | `runs/fish_te_axis_treatment_20261001/tables/06_provenance.json` | differ | text cells differing: 6 | run metadata, as above |
| `2026-10-01_te_axis_treatment/tables/06_te_axis_qc_by_treatment.tsv` | `runs/fish_te_axis_treatment_20261001/tables/06_te_axis_qc_by_treatment.tsv` | differ | 164 of 216 cells differ; max abs 4.22e-15, max rel 3.06e-13 |  |
| `2026-10-01_te_axis_treatment/tables/06_te_axis_scores.tsv` | `runs/fish_te_axis_treatment_20261001/tables/06_te_axis_scores.tsv` | differ | 146 of 1296 cells differ; max abs 9.77e-15, max rel 3.54e-13 |  |
| `2026-10-01_te_axis_treatment/tables/06_te_axis_treatment_no_outliers.tsv` | none | not produced | | not made by `compute` |
| `2026-10-01_te_axis_treatment/tables/06_te_axis_treatment_tests.tsv` | `runs/fish_te_axis_treatment_20261001/tables/06_te_axis_treatment_tests.tsv` | differ | 810 of 1134 cells differ; max abs 1.65e-13, max rel 9.2e-12 |  |
| `2026-10-01_te_axis_treatment/tables/06_te_read_pct.tsv` | `runs/fish_te_axis_treatment_20261001/tables/06_te_read_pct.tsv` | differ | 75 of 108 cells differ; max abs 2.7e-13, max rel 1.56e-13 |  |
| `2026-10-01_zebrafish_validation/tables/08_class_te_response.tsv` | `runs/zf_validation_20261001/tables/08_class_te_response.tsv` | differ | 55 of 309 cells differ; max abs 1.11e-16, max rel 8.5e-16 |  |
| `2026-10-01_zebrafish_validation/tables/08_family_lfc.tsv.gz` | `runs/zf_validation_20261001/tables/08_family_lfc.tsv.gz` | differ | 0 of 46173 cells differ; max abs 0, max rel 0 | content identical (0 of 46,173 cells differ); the gzip header differs |
| `2026-10-01_zebrafish_validation/tables/08_perturbed_gene_check.tsv` | `runs/zf_validation_20261001/tables/08_perturbed_gene_check.tsv` | identical |  |  |
| `2026-10-01_zebrafish_validation/tables/08_prediction_concordance.tsv` | `runs/zf_validation_20261001/tables/08_prediction_concordance.tsv` | identical |  |  |
| `2026-10-01_zebrafish_validation/tables/08_study_te_response.tsv` | `runs/zf_validation_20261001/tables/08_study_te_response.tsv` | differ | 35 of 490 cells differ; max abs 3.41e-13, max rel 2.89e-15 |  |
| `2026-10-01_zebrafish_validation/tables/08c_qc_by_arm.tsv` | `runs/zf_validation_20261001/tables/08c_qc_by_arm.tsv` | differ | 73 of 442 cells differ; max abs 1.14e-13, max rel 4.39e-14 |  |
| `2026-10-01_zebrafish_validation/tables/08c_qc_per_run.tsv` | `runs/zf_validation_20261001/tables/08c_qc_per_run.tsv` | identical |  |  |
| `2026-10-01_zebrafish_validation/tables/08d_hUHRF1_OE_te_share_contributors.tsv` | `runs/zf_validation_20261001/tables/08d_hUHRF1_OE_te_share_contributors.tsv` | differ | 0 of 4080 cells differ; max abs 0, max rel 0 | same values; the order of tied rows (value 0.0) differs |
| `2026-10-01_zebrafish_validation/tables/08d_hUHRF1_OE_te_share_contributors_read_totals.tsv` | `runs/zf_validation_20261001/tables/08d_hUHRF1_OE_te_share_contributors_read_totals.tsv` | differ | 7 of 30 cells differ; max abs 1.78e-15, max rel 1.25e-15 |  |
| `2026-10-01_zebrafish_validation/tables/08e_clutch_deviation.tsv` | `runs/zf_validation_20261001/tables/08e_clutch_deviation.tsv` | identical |  | not recomputed: the committed table was copied in as the input of 08f |
| `2026-10-01_zebrafish_validation/tables/08f_exploratory_concordance.tsv` | `runs/zf_validation_20261001/tables/08f_exploratory_concordance.tsv` | identical |  | recomputed from the copied 08e table and the rerun 08 tables |
| `2026-10-01_zebrafish_validation/tables/08s_gse234993_mapping.tsv` | none | not produced | | not made by `compute` |
| `2026-10-01_zebrafish_validation/tables/08s_gse234993_mapping_tests.tsv` | none | not produced | | not made by `compute` |
| `2026-10-01_zebrafish_validation/tables/08s_gse234993_sensitivity.tsv` | none | not produced | | not made by `compute` |

## Limitations

- Byte identity was not expected for tables that print floating-point values at full precision, and we did not get it. The numeric comparison is the useful one for those files; it reads each table with pandas and compares cell by cell.
- The fixed `run_all.sh` was not run end to end. Its commands equal the ones we ran in steps 1 to 4 above.
- 08e is not recomputed, so 08f checks only its own code on the committed 08e table.
- One machine, one environment. We did not repeat the run, so we cannot say whether the boundary behaviour in `06_phenotype_tests.tsv` also changes between runs on the same machine.
