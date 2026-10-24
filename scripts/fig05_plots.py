"""Figure 05_fig1_real_vs_random: the fish bulk ISP regulators against expression-matched random genes.

Added in this repository (the 2026-10-01 run drew no figure). It plots only values saved in
05c_*_per_regulator.tsv and 05c_*_overall.tsv; no number is recomputed. Delete arm.
"""
from pathlib import Path

import numpy as np
import pandas as pd


def plot_05(tables: Path, figs: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcdefaults()
    plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})
    modes = [("plain", "Plain (within tissue x treatment)"), ("rmTEpc1", "Global TE axis removed")]
    fig, axes = plt.subplots(2, 2, figsize=(8.6, 7.2))
    for j, (mode, title) in enumerate(modes):
        per = pd.read_csv(tables / f"05c_{mode}_per_regulator.tsv", sep="\t")
        ov = pd.read_csv(tables / f"05c_{mode}_overall.tsv", sep="\t")
        per = per[per.perturb_type == "delete"].sort_values("p_emp_hits", ascending=False).reset_index(drop=True)
        ov = ov[ov.perturb_type == "delete"].iloc[0]
        y = np.arange(len(per))
        a = axes[0, j]
        a.scatter(per.n_fdr05, y, s=14, color="#2166ac")
        a.axvline(ov.median_hits_random, color="0.4", ls="--", lw=0.8)
        a.axvline(ov.median_hits_real, color="#2166ac", ls=":", lw=0.8)
        a.set_yticks(y, per.target, fontsize=6.5)
        a.set_xlabel("TE families with FDR < 0.05 (of those tested)")
        a.set_title(f"{title}\nmedian: real {ov.median_hits_real:g} (dotted), random {ov.median_hits_random:g} (dashed); "
                    f"Mann-Whitney p = {ov.mw_p_n_fdr05:.2g}", fontsize=7.5)
        b = axes[1, j]
        b.scatter(per.p_emp_hits, y, s=14, color="#b2182b")
        b.axvline(0.05, color="0.4", ls="--", lw=0.8)
        b.set_xlim(0, 1.02)
        b.set_yticks(y, per.target, fontsize=6.5)
        b.set_xlabel("empirical p against random pseudo-regulators (hits)")
        b.set_title(f"{int(ov['n_regulators_p_emp_hits_le_0.05'])} of {int(ov.n_regulators)} regulators at p <= 0.05 "
                    f"(expected {ov.expected_by_chance:g}); {int(ov['n_regulators_fdr_hits_lt_0.05'])} pass BH", fontsize=7.5)
    fig.suptitle("Fish bulk ISP, delete arm: 30 TE regulators against 90 expression-matched random genes "
                 "(A. polyacanthus, 162 libraries)", fontsize=8.5)
    fig.tight_layout()
    fig.savefig(figs / "05_fig1_real_vs_random.png", dpi=300)
    fig.savefig(figs / "05_fig1_real_vs_random.pdf")
    plt.close(fig)
    plt.rcdefaults()  # do not leak these settings into the figures drawn after this one
