#!/usr/bin/env python3
"""05c — Are real regulators more TE-informative than random pseudo-regulators?

Compares a 05_bulk_linear_isp.py run on the configured regulators with one
or more `--random-regulators` runs (expression-matched random genes), per
perturbation arm, with statistics:

  per regulator
    n_fdr05             TE families with FDR < 0.05
    median_abs_effect   median |effect_vs_null| over tested families
    p_emp_hits          (1 + #random >= n_fdr05) / (1 + #random), one-sided
    p_emp_effect        the same for median_abs_effect
    fdr_hits/effect     BH across regulators
    frac_up_sig         share of significant families predicted to go up,
                        with an exact binomial 95% CI
  overall
    Mann-Whitney U (two-sided): real vs random n_fdr05 and median_abs_effect
    binomial test: number of regulators with p_emp_hits <= 0.05 vs the 5%
    expected by chance

Usage:
  python scripts/05c_compare_random_controls.py --real <05_..._effects.tsv> \
      --random <05_random1_..._effects.tsv> [...] --out <dir> [--label plain]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

FDR = 0.05


def per_regulator(df: pd.DataFrame) -> pd.DataFrame:
    t = df[df["status"] == "tested"]
    rows = []
    for (reg, arm), g in t.groupby(["target", "perturb_type"]):
        sig = g[g["fdr"] < FDR]
        n_up = int((sig["effect_vs_null"] > 0).sum())
        ci = (stats.binomtest(n_up, len(sig)).proportion_ci() if len(sig) else None)
        rows.append({"target": reg, "perturb_type": arm, "n_tested": len(g), "n_fdr05": len(sig),
                     "median_abs_effect": float(np.median(np.abs(g["effect_vs_null"]))),
                     "frac_up_sig": n_up / len(sig) if len(sig) else np.nan,
                     "frac_up_ci95_low": ci.low if ci else np.nan,
                     "frac_up_ci95_high": ci.high if ci else np.nan})
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--real", required=True)
    ap.add_argument("--random", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--label", default="run")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    real = per_regulator(pd.read_csv(args.real, sep="\t"))
    rand = pd.concat([per_regulator(pd.read_csv(f, sep="\t")).assign(seed=Path(f).name)
                      for f in args.random], ignore_index=True)
    per, overall = [], []
    for arm, r in real.groupby("perturb_type"):
        q = rand[rand["perturb_type"] == arm]
        r = r.copy()
        for col, name in (("n_fdr05", "hits"), ("median_abs_effect", "effect")):
            r[f"p_emp_{name}"] = [(1 + (q[col] >= v).sum()) / (1 + len(q)) for v in r[col]]
            r[f"fdr_{name}"] = multipletests(r[f"p_emp_{name}"], method="fdr_bh")[1]
        per.append(r)
        k = int((r["p_emp_hits"] <= 0.05).sum())
        rec = {"label": args.label, "perturb_type": arm, "n_regulators": len(r), "n_random": len(q),
               "median_hits_real": float(r["n_fdr05"].median()),
               "median_hits_random": float(q["n_fdr05"].median()),
               "median_abs_effect_real": float(r["median_abs_effect"].median()),
               "median_abs_effect_random": float(q["median_abs_effect"].median()),
               "n_regulators_p_emp_hits_le_0.05": k,
               "expected_by_chance": 0.05 * len(r),
               "p_binom_excess": float(stats.binomtest(k, len(r), 0.05, alternative="greater").pvalue),
               "n_regulators_fdr_hits_lt_0.05": int((r["fdr_hits"] < 0.05).sum())}
        for col in ("n_fdr05", "median_abs_effect"):
            mw = stats.mannwhitneyu(r[col], q[col], alternative="two-sided")
            rec[f"mw_U_{col}"], rec[f"mw_p_{col}"] = float(mw.statistic), float(mw.pvalue)
        overall.append(rec)
    pd.concat(per).sort_values(["perturb_type", "p_emp_hits"]).to_csv(
        out / f"05c_{args.label}_per_regulator.tsv", sep="\t", index=False)
    ov = pd.DataFrame(overall)
    ov.to_csv(out / f"05c_{args.label}_overall.tsv", sep="\t", index=False)
    print(ov.T.to_string())


if __name__ == "__main__":
    main()
