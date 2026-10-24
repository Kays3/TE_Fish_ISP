"""Summarize fish bulk ISP: real regulators vs expression-matched random genes."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

root = Path(sys.argv[1])            # local folder holding pulled tables
mode = sys.argv[2] if len(sys.argv) > 2 else ""     # "" or "rmTEpc1_"
real = pd.read_csv(root / f"real/05_{mode}bulk_isp_effects.tsv", sep="\t")
pattern = f"random*/05_random[0-9]_{mode}bulk_isp_effects.tsv"
rand = pd.concat([pd.read_csv(f, sep="\t").assign(seed=f.parent.name)
                  for f in sorted(root.glob(pattern))], ignore_index=True)
print(f"mode={mode or 'plain'}; random seeds: {rand.seed.nunique()}, pseudo-regulators: {rand.target.nunique()}")
FDR = 0.05


def per_regulator(df):
    t = df[df["status"] == "tested"]
    g = t.groupby(["target", "perturb_type"])
    return pd.DataFrame({
        "n_tested": g.size(),
        "n_fdr05": g["fdr"].apply(lambda x: int((x < FDR).sum())),
        "median_abs_effect": g["effect_vs_null"].apply(lambda x: float(np.median(np.abs(x)))),
        "frac_up_among_sig": g[['fdr','effect_vs_null']].apply(lambda d: float((d.loc[d.fdr < FDR, "effect_vs_null"] > 0).mean())
                                     if (d.fdr < FDR).any() else np.nan),
    }).reset_index()


pr = per_regulator(real)
rr = per_regulator(rand)
print("== real regulators, delete arm ==")
print(pr[pr.perturb_type == "delete"].sort_values("n_fdr05", ascending=False).to_string(index=False))
print("\n== random pseudo-regulators, delete arm (all seeds) ==")
print(rr[rr.perturb_type == "delete"].describe().round(3).to_string())

d_real = pr[pr.perturb_type == "delete"]
d_rand = rr[rr.perturb_type == "delete"]
print("\nmedian per-regulator n_fdr05: real %.1f vs random %.1f" % (d_real.n_fdr05.median(), d_rand.n_fdr05.median()))
print("median per-regulator median|effect|: real %.4f vs random %.4f" % (
    d_real.median_abs_effect.median(), d_rand.median_abs_effect.median()))
# empirical: for each real regulator, fraction of random pseudo-regulators with >= as many hits
d_real = d_real.assign(
    pct_random_ge_hits=[float((d_rand.n_fdr05 >= n).mean()) for n in d_real.n_fdr05],
    pct_random_ge_effect=[float((d_rand.median_abs_effect >= e).mean()) for e in d_real.median_abs_effect])
print("\n== real vs random calibration (delete) ==")
print(d_real.sort_values("pct_random_ge_hits")[["target", "n_fdr05", "median_abs_effect",
      "frac_up_among_sig", "pct_random_ge_hits", "pct_random_ge_effect"]].round(3).to_string(index=False))
d_real.to_csv(root / f"real_vs_random_delete_{mode or 'plain_'}.tsv", sep="\t", index=False)

# class-level view for the regulators that beat the random control
keep = d_real[d_real.pct_random_ge_hits <= 0.05].target
sig = real[(real.status == "tested") & (real.fdr < FDR) & real.target.isin(keep) & (real.perturb_type == "delete")]
if len(sig):
    cls = sig.assign(up=sig.effect_vs_null > 0).groupby(["target", "te_class"]).agg(
        n=("up", "size"), frac_up=("up", "mean")).reset_index()
    print("\n== TE classes hit by regulators beating the control (delete) ==")
    print(cls.sort_values(["target", "n"], ascending=[True, False]).round(2).to_string(index=False))
    cls.to_csv(root / f"class_summary_delete_{mode or 'plain_'}.tsv", sep="\t", index=False)
