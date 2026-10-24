# Figures index

Ten figures are committed, each as `.png` and `.pdf`, byte-for-byte as in `Kays3/Geneformer_TE` at `0f0d812`
(`docs/runs/<run>/figures/`). `scripts/make_figures.py` rebuilds them from the committed tables only; no counts, no
NAS and no GPU are needed:

```bash
PYTHON=<pinned env python> scripts/run_all.sh figures      # = make_figures.py --out rebuilt --check
```

`--check` compares each rebuilt PNG with the committed one by sha256. With the pinned environment
(`environment/README.md`, matplotlib 3.11.2) all ten are identical; with another matplotlib minor version the data are
the same and the layout shifts. PDFs are rebuilt too but carry a creation date, so they are not compared by bytes.
Checksums of the committed PNGs are in `manifests/figures.sha256`.

The 2026-10-01 fish bulk ISP run (05) drew no figure, and none is added here: this repository makes no new analysis.

| # | Figure (committed PNG and PDF) | Drawn from the tables of its run | Code | PNG sha256 (first 16) |
|---|---|---|---|---|
| 1 | [`06_te_axis_by_treatment`](../docs/runs/2026-10-01_te_axis_treatment/figures/06_te_axis_by_treatment.png): Global TE axis and strict-TE read fraction by CO2 treatment, per tissue, with family-adjusted contrasts and BH-FDR | `06_te_axis_scores.tsv, 06_te_axis_treatment_tests.tsv` | `scripts/fig06_plots.py` | `b031fe95813ae809` |
| 2 | [`06_phenotype_contrasts`](../docs/runs/2026-10-01_te_axis_treatment/figures/06_phenotype_contrasts.png): Tolerant minus sensitive parental phenotype within tissue (mixed model), three TE measures | `06_phenotype_tests.tsv` | `scripts/fig06_plots.py` | `dba2e2a4181a0a73` |
| 3 | [`07_fig1_gene_effects_heatmap`](../docs/runs/2026-10-01_gene_axis/figures/07_fig1_gene_effects_heatmap.png): Treatment effect on each candidate gene (log2FC) by tissue and contrast | `07_gene_treatment_effects.tsv` | `scripts/07_treatment_gene_axis.py, make_figures` | `7496efdfc3fadcbf` |
| 4 | [`07_fig2_set_shift`](../docs/runs/2026-10-01_gene_axis/figures/07_fig2_set_shift.png): Set-level treatment shift against expression-matched random gene sets | `07_set_treatment_shift.tsv` | `scripts/07_treatment_gene_axis.py, make_figures` | `91f33453b078f831` |
| 5 | [`07_fig3_gene_mediation`](../docs/runs/2026-10-01_gene_axis/figures/07_fig3_gene_mediation.png): Gene-level mediation of the treatment effect on the TE axis | `07_gene_mediation.tsv` | `scripts/07_treatment_gene_axis.py, make_figures` | `11d6e7b77a75d5a4` |
| 6 | [`07_fig4_axis_association`](../docs/runs/2026-10-01_gene_axis/figures/07_fig4_axis_association.png): Within-group association of each gene with the TE axis | `07_gene_axis_association.tsv` | `scripts/07_treatment_gene_axis.py, make_figures` | `06762b0bb343e192` |
| 7 | [`07_fig5_set_mediation`](../docs/runs/2026-10-01_gene_axis/figures/07_fig5_set_mediation.png): Set-level mediation against random sets of the same size | `07_set_mediation.tsv` | `scripts/07_treatment_gene_axis.py, make_figures` | `ebbbd70648b58775` |
| 8 | [`08_fig1_study_te_response`](../docs/runs/2026-10-01_zebrafish_validation/figures/08_fig1_study_te_response.png): Zebrafish: change in TE share per study with 95% CI, and TE families called up and down | `08_study_te_response.tsv` | `scripts/08_zebrafish_validation.py, figures` | `c1751266c594f335` |
| 9 | [`08_fig2_class_response`](../docs/runs/2026-10-01_zebrafish_validation/figures/08_fig2_class_response.png): Zebrafish: median log2FC of TE families per class and study | `08_class_te_response.tsv` | `scripts/08_zebrafish_validation.py, figures` | `21ff190c84dca118` |
| 10 | [`08_fig3_prediction_vs_observation`](../docs/runs/2026-10-01_zebrafish_validation/figures/08_fig3_prediction_vs_observation.png): Fish ISP prediction against the zebrafish TE-share change (Spearman rho, n = 17) | `08_study_te_response.tsv (fish_*_pred_frac_up columns)` | `scripts/08_zebrafish_validation.py, figures` | `1ff9f5cb3729d95f` |

## Figures whose inputs are missing

None. All ten rebuild from the committed tables.

## Reading the figures (ISP-STD-1 wording)

- Figures 8 to 10: the zebrafish studies are the unit. A study with a dagger has a genotype-linked mapping gap.
  Asterisks mark BH-FDR < 0.05 across studies; no loss-of-function study shifts global TE share at that threshold.
- Figure 10 crosses species, tissue and perturbation type, and several zebrafish studies share one fish prediction.
- Figures 3 to 7 are exploratory candidate-gene analyses judged against expression-matched random gene sets.
