# Figures index

Every figure is rebuilt by `scripts/make_figures.py` (or `scripts/run_all.sh figures`) from the saved tables in
`results/`; no counts and no GPU are needed. Each figure has a `.png` and a `.pdf` next to it. Checksums of the PNGs are
in `manifests/figures.sha256`; PDFs carry a fixed creation date and also repeat byte for byte.

The committed files were made with the pinned environment (`environment/README.md`: matplotlib 3.11.2). Figures 2
to 11 are byte-identical to the PNGs committed in `Kays3/Geneformer_TE` (`docs/runs/2026-10-01_*/figures/`),
checked on 25 October 2026 against `origin/main` at `0f0d812`. Figure 1 did not exist there; it plots values saved in
the 05c tables only.

| # | Figure | Drawn from `results/<run>/tables/` | Code | PNG sha256 (first 16) | Source PNG in Geneformer_TE |
|---|---|---|---|---|---|
| 1 | [`05_fig1_real_vs_random`](2026-10-01_fish_bulk_isp/05_fig1_real_vs_random.png) (2026-10-01_fish_bulk_isp): Fish bulk ISP, delete arm: each of the 30 regulators against 90 expression-matched random genes, plain and with the TE axis removed (families at FDR < 0.05; empirical p) | 05c_plain_per_regulator.tsv, 05c_rmTEpc1_per_regulator.tsv, 05c_plain_overall.tsv, 05c_rmTEpc1_overall.tsv | scripts/fig05_plots.py | `7f8f43dcddb1c5e5` | none: new here |
| 2 | [`06_te_axis_by_treatment`](2026-10-01_te_axis_treatment/06_te_axis_by_treatment.png) (2026-10-01_te_axis_treatment): Global TE axis and strict-TE read fraction by CO2 treatment, per tissue, with family-adjusted contrasts and BH-FDR | 06_te_axis_scores.tsv, 06_te_axis_treatment_tests.tsv | scripts/fig06_plots.py | `b031fe95813ae809` | byte-identical |
| 3 | [`06_phenotype_contrasts`](2026-10-01_te_axis_treatment/06_phenotype_contrasts.png) (2026-10-01_te_axis_treatment): Tolerant minus sensitive parental phenotype within tissue (mixed model), three TE measures | 06_phenotype_tests.tsv | scripts/fig06_plots.py | `dba2e2a4181a0a73` | byte-identical |
| 4 | [`07_fig1_gene_effects_heatmap`](2026-10-01_gene_axis/07_fig1_gene_effects_heatmap.png) (2026-10-01_gene_axis): Treatment effect on each candidate gene (log2FC) by tissue and contrast | 07_gene_treatment_effects.tsv | scripts/07_treatment_gene_axis.py, make_figures | `7496efdfc3fadcbf` | byte-identical |
| 5 | [`07_fig2_set_shift`](2026-10-01_gene_axis/07_fig2_set_shift.png) (2026-10-01_gene_axis): Set-level treatment shift against expression-matched random gene sets | 07_set_treatment_shift.tsv | scripts/07_treatment_gene_axis.py, make_figures | `91f33453b078f831` | byte-identical |
| 6 | [`07_fig3_gene_mediation`](2026-10-01_gene_axis/07_fig3_gene_mediation.png) (2026-10-01_gene_axis): Gene-level mediation of the treatment effect on the TE axis | 07_gene_mediation.tsv | scripts/07_treatment_gene_axis.py, make_figures | `11d6e7b77a75d5a4` | byte-identical |
| 7 | [`07_fig4_axis_association`](2026-10-01_gene_axis/07_fig4_axis_association.png) (2026-10-01_gene_axis): Within-group association of each gene with the TE axis | 07_gene_axis_association.tsv | scripts/07_treatment_gene_axis.py, make_figures | `06762b0bb343e192` | byte-identical |
| 8 | [`07_fig5_set_mediation`](2026-10-01_gene_axis/07_fig5_set_mediation.png) (2026-10-01_gene_axis): Set-level mediation against random sets of the same size | 07_set_mediation.tsv | scripts/07_treatment_gene_axis.py, make_figures | `ebbbd70648b58775` | byte-identical |
| 9 | [`08_fig1_study_te_response`](2026-10-01_zebrafish_validation/08_fig1_study_te_response.png) (2026-10-01_zebrafish_validation): Zebrafish: change in TE share per study with 95% CI, and TE families called up and down | 08_study_te_response.tsv | scripts/08_zebrafish_validation.py, figures | `c1751266c594f335` | byte-identical |
| 10 | [`08_fig2_class_response`](2026-10-01_zebrafish_validation/08_fig2_class_response.png) (2026-10-01_zebrafish_validation): Zebrafish: median log2FC of TE families per class and study | 08_class_te_response.tsv | scripts/08_zebrafish_validation.py, figures | `21ff190c84dca118` | byte-identical |
| 11 | [`08_fig3_prediction_vs_observation`](2026-10-01_zebrafish_validation/08_fig3_prediction_vs_observation.png) (2026-10-01_zebrafish_validation): Fish ISP prediction against the zebrafish TE-share change (Spearman rho, n = 17) | 08_study_te_response.tsv (fish_*_pred_frac_up columns) | scripts/08_zebrafish_validation.py, figures | `1ff9f5cb3729d95f` | byte-identical |

## Figures whose inputs are missing

None. All eleven rebuild from the committed tables. Two figures of the source work are not rebuilt because they
never existed as files: the family-level concordance and timeline panels of the history page were drawn inside that
HTML page (the concordance values are in `results/2026-10-01_zebrafish_validation/tables/08f_exploratory_concordance.tsv`).

## Reading the figures (ISP-STD-1 wording)

- Figure 1: the regulators are compared with expression-matched random genes (B1). One regulator reaches an
  empirical p <= 0.05 in the plain run (1.5 expected by chance) and none in the axis-removed run; none passes BH.
- Figures 9 to 11: the zebrafish studies are the unit. A study with a dagger has a genotype-linked mapping gap.
  Asterisks mark BH-FDR < 0.05 across studies. No loss-of-function study shifts global TE share at that threshold.
- Figure 11 crosses species, tissue and perturbation type, and several zebrafish studies share one fish prediction.
