#!/usr/bin/env python3
"""Exploratory cross-study concordance on the family-level lean (not prespecified).

Unit: a study whose significant TE families lean one way (binomial FDR < 0.05 in
08_study_te_response.tsv). Rows: silencer loss as configured (every KO/KD), the
DNA-methylation subset (dnmt1, uhrf1, dnmt3), that subset without QC-flagged
studies (08c flag), and the two fish bulk-ISP directions. "planned" replaces the
three clutch-paired studies with their ~ role result (08e). Exact two-sided
binomial vs 0.5 with Clopper-Pearson 95% CI; BH-FDR across rows. Also the two
Spearman correlations of fish prediction with TE-share change.

Usage: 08f_exploratory_concordance.py TABLES_DIR
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

DNA_METH = {"DNMT1", "UHRF1", "DNMT3B", "DNMT3A"}


def main(t: Path) -> None:
    S = pd.read_csv(t / "08_study_te_response.tsv", sep="\t")
    Q = pd.read_csv(t / "08c_qc_by_arm.tsv", sep="\t").set_index("study_id")
    E = pd.read_csv(t / "08e_clutch_deviation.tsv", sep="\t")
    S["lean"] = np.where(S.fdr_binom < 0.05, np.sign(S.frac_up_sig - 0.5), 0)
    S["flag"] = S.study_id.map(Q["flag_mapping_gap"]).fillna(False).astype(bool)
    planned = S.copy()
    for r in E[E.design == "~ role"].itertuples():
        lean = np.sign(r.frac_up - 0.5) if r.fdr_bh_binom < 0.05 and np.isfinite(r.frac_up) else 0
        planned.loc[planned.study_id == r.study_id, "lean"] = lean
    lof = lambda d: d[d.perturbation.isin(["KO", "KD"]) & (d.lean != 0)]  # noqa: E731
    rows = []

    def add(label, design, d, pred):
        k = int((d.lean == pred(d)).sum())
        n = len(d)
        bt = stats.binomtest(k, n, 0.5) if n else None
        ci = bt.proportion_ci(method="exact") if bt else None
        rows.append({"prediction": label, "design": design, "n_studies": n, "n_agree": k,
                     "frac_agree": k / n if n else np.nan, "ci95_low": ci.low if ci else np.nan,
                     "ci95_high": ci.high if ci else np.nan, "p_binom": bt.pvalue if bt else np.nan,
                     "n_gse": d.gse.nunique()})

    up = lambda d: 1  # noqa: E731
    for design, D in (("primary", S), ("planned", planned)):
        add("silencer loss -> TE up (all KO/KD)", design, lof(D), up)
        dm = lof(D)[lof(D).target.isin(DNA_METH)]
        add("DNA-methylation loss -> TE up", design, dm, up)
        add("DNA-methylation loss -> TE up, no QC flag", design, dm[~dm.flag], up)
    lean = S[S.lean != 0]
    for m in ("plain", "rmTEpc1"):
        d = lean[lean[f"fish_{m}_pred_direction"] != 0]
        add(f"fish bulk ISP ({m})", "primary", d, lambda x, m=m: x[f"fish_{m}_pred_direction"])
    T = pd.DataFrame(rows)
    T["fdr_bh"] = np.nan
    ok = T.p_binom.notna()
    T.loc[ok, "fdr_bh"] = multipletests(T.loc[ok, "p_binom"], method="fdr_bh")[1]
    for m in ("plain", "rmTEpc1"):
        d = S.dropna(subset=[f"fish_{m}_pred_frac_up"])
        rho, p = stats.spearmanr(d[f"fish_{m}_pred_frac_up"], d.te_share_rel_change_pct)
        T = pd.concat([T, pd.DataFrame([{"prediction": f"Spearman fish ({m}) vs TE-share change",
                                          "design": "primary", "n_studies": len(d), "rho": rho, "p_spearman": p}])])
    T.to_csv(t / "08f_exploratory_concordance.tsv", sep="\t", index=False)
    print(T.round(3).to_string(index=False))


if __name__ == "__main__":
    main(Path(sys.argv[1]))
