"""Figures 06_te_axis_by_treatment and 06_phenotype_contrasts, rebuilt from the saved 06 tables.

The plotting code is the block labelled "# Figure" in 06_te_axis_treatment.py (copied without
change; that script draws the figures inline, so they could not be rebuilt without the counts).
Inputs: 06_te_axis_scores.tsv (sc), 06_te_axis_treatment_tests.tsv (res), 06_phenotype_tests.tsv (pres).
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


def plot_06(sc: pd.DataFrame, res: pd.DataFrame, pres: pd.DataFrame, figs: Path) -> None:
    # Figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    order = ["control", "developmental", "intergenerational"]
    fig, axes = plt.subplots(2, 3, figsize=(11, 6.5), sharey="row")
    rng = np.random.default_rng(0)
    for r, measure in enumerate(["axis_score", "te_frac_strict"]):
        for c, tissue in enumerate(["brain", "gill", "liver"]):
            ax = axes[r, c]
            d = sc[sc["tissue"] == tissue]
            for k, (pheno, col, off) in enumerate([("sensitive", "#c0504d", -0.15), ("tolerant", "#4f81bd", 0.15)]):
                for i, t in enumerate(order):
                    v = d[(d["treatment"] == t) & (d["parental_phenotype"] == pheno)][measure]
                    ax.scatter(i + off + rng.uniform(-0.05, 0.05, len(v)), v, s=12, color=col, alpha=0.7,
                               label=pheno if (i == 0 and c == 0 and r == 0) else None)
                    ci = float(stats.t.ppf(0.975, len(v) - 1)) * v.std(ddof=1) / np.sqrt(len(v))
                    ax.errorbar(i + off, v.mean(), yerr=ci, color=col, lw=1.5, capsize=3, zorder=3)
                    ax.hlines(v.mean(), i + off - 0.1, i + off + 0.1, color=col, lw=2)
            # model-based contrasts (primary pooled model, family-adjusted) with BH-FDR
            tm = "axis_score" if measure == "axis_score" else "te_frac_strict"
            lines = []
            for cc in res[(res["model"] == "tissue_pooled") & (res["measure"] == tm)
                          & (res["tissue"] == tissue)].itertuples():
                lab = "dev" if cc.contrast.startswith("developmental") else "intergen"
                star = "*" if cc.fdr < 0.05 else ""
                lines.append(f"{lab}: {cc.estimate:+.2f} [{cc.ci95_low:+.2f}, {cc.ci95_high:+.2f}] "
                             f"FDR {cc.fdr:.2g}{star}")
            ax.text(0.02, 0.02, "\n".join(lines), transform=ax.transAxes, fontsize=6.5,
                    va="bottom", ha="left", family="monospace",
                    bbox=dict(facecolor="white", alpha=0.8, edgecolor="none", pad=1.5))
            ax.set_xticks(range(3), ["ctrl", "dev", "intergen"])
            if r == 0:
                ax.set_title(tissue)
            if c == 0:
                ax.set_ylabel({"axis_score": "Global TE axis (SD)",
                               "te_frac_strict": "logit TE read fraction"}[measure])
    for row in axes:                          # room below the data for the statistics text
        lo, hi = row[0].get_ylim()
        row[0].set_ylim(lo - 0.32 * (hi - lo), hi)
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.suptitle("Global TE-load axis by CO2 treatment (A. polyacanthus, 162 libraries)\n"
                 "points = libraries; bars = group mean ± 95% CI; text = family-adjusted contrast vs "
                 "control [95% CI], BH-FDR (* < 0.05)", fontsize=9)
    fig.tight_layout()
    fig.savefig(figs / "06_te_axis_by_treatment.png", dpi=200)
    fig.savefig(figs / "06_te_axis_by_treatment.pdf")
    # Forest plot: tolerant - sensitive (mixed model) per tissue x treatment
    fig2, axes2 = plt.subplots(1, 3, figsize=(14, 5.2), sharey=True)
    labels = {"tolerant_minus_sensitive@control": "ctrl", "tolerant_minus_sensitive@developmental": "dev",
              "tolerant_minus_sensitive@intergenerational": "intergen",
              "tolerant_minus_sensitive@mean_of_treatments": "mean",
              "interaction_developmental_effect_tol_minus_sens": "Δdev effect",
              "interaction_intergenerational_effect_tol_minus_sens": "Δintergen effect"}
    ylabels = []
    for panel, (ax, measure) in enumerate(zip(axes2, ("axis_score", "te_frac_strict", "te_frac_all"))):
        sub = pres[(pres["method"] == "mixed") & (pres["measure"] == measure)]
        for k, (tissue, col) in enumerate([("brain", "#7f6aa8"), ("gill", "#2a9d8f"), ("liver", "#c97b2a")]):
            for j, (cname, lab) in enumerate(labels.items()):
                r = sub[(sub["tissue"] == tissue) & (sub["contrast"] == cname)].iloc[0]
                y = k * (len(labels) + 1) + j
                ax.errorbar(r.estimate, y, xerr=[[r.estimate - r.ci95_low], [r.ci95_high - r.estimate]],
                            fmt="o", color=col, ms=4, capsize=2)
                ax.text(1.02, y, f"p={r.p:.2g}" + ("*" if r.fdr < 0.05 else ""),
                        transform=ax.get_yaxis_transform(), fontsize=6.5, va="center", ha="left",
                        clip_on=False)
                if panel == 0:
                    ylabels.append((y, f"{tissue} {lab}"))
        ax.axvline(0, color="k", lw=0.7)
        ax.set_title({"axis_score": "Global TE axis (SD)", "te_frac_strict": "logit strict-TE fraction",
                      "te_frac_all": "logit all-repeat fraction"}[measure], fontsize=9)
    axes2[0].set_yticks([y for y, _ in ylabels])
    axes2[0].set_yticklabels([t for _, t in ylabels], fontsize=7)
    axes2[0].tick_params(axis="y", labelleft=True)
    axes2[0].set_ylim(max(y for y, _ in ylabels) + 0.8, -0.8)   # top-to-bottom order
    fig2.suptitle("Tolerant − sensitive parental phenotype within tissue: estimate ± 95% CI from a mixed model "
                  "with family random intercept\n(fish-level OLS where the family variance is at its boundary); "
                  "raw p shown, * BH-FDR < 0.05. Δ = phenotype × treatment interaction", fontsize=9)
    fig2.tight_layout(rect=(0, 0, 0.97, 1))
    fig2.subplots_adjust(wspace=0.32)
    fig2.savefig(figs / "06_phenotype_contrasts.png", dpi=200)
    fig2.savefig(figs / "06_phenotype_contrasts.pdf")
