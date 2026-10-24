#!/usr/bin/env python3
"""Which TE families carry a study's change in TE share of counted reads
(review item: hUHRF1 OE share up while families lean down).

For each strict TE family: mean share of counted reads per arm, the
perturbed - control difference, and its fraction of the total TE-share change.
Also reports the arm means of TE and gene counted reads (millions), so a
composition effect (genes falling) can be told from TE reads rising.

Usage: 08d_te_share_contributors.py COUNTS_DIR SAMPLES_TSV STUDY_ID OUT_TSV [QC_DIR]
With QC_DIR (STAR Log.final.out files), also writes <OUT_TSV stem>_read_totals.tsv:
input, uniquely mapped, TE- and gene-assigned reads per arm with Welch p.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils.common import NON_TE_CLASSES  # noqa: E402


def main(counts_dir: Path, samples: Path, study: str, out: Path, qc_dir: Path | None = None) -> None:
    S = pd.read_csv(samples, sep="\t")
    S = S[S.parental_phenotype == study]
    M = pd.concat({r.sample_id: pd.read_csv(counts_dir / f"{r.sample_id}.cntTable", sep="\t", index_col=0).iloc[:, 0]
                   for r in S.itertuples()}, axis=1).fillna(0)
    parts = M.index.to_series().str.split(":")
    cls = parts.str[2].where(parts.str.len() >= 3).str.rstrip("?")
    strict = cls.notna() & ~cls.isin(NON_TE_CLASSES | {"Satellite", "Low_complexity"})
    tot = M.sum()
    share = M[strict] / tot
    pert = S.loc[S.treatment == "perturbed", "sample_id"]
    ctrl = S.loc[S.treatment == "control", "sample_id"]
    d = share[pert].mean(axis=1) - share[ctrl].mean(axis=1)
    total = d.sum()
    F = pd.DataFrame({"family": M.index[strict], "te_class": cls[strict].values,
                      "share_control_pct": 100 * share[ctrl].mean(axis=1).values,
                      "share_perturbed_pct": 100 * share[pert].mean(axis=1).values,
                      "delta_share_pct_points": 100 * d.values,
                      "frac_of_total_delta": (d / total).values})
    F = F.sort_values("delta_share_pct_points", ascending=False)
    F.to_csv(out, sep="\t", index=False)
    te_m, gene_m = M[strict].sum() / 1e6, M[~cls.notna()].sum() / 1e6
    print(f"{study}: TE share change {100 * total:+.3f} points")
    series = [("TE-assigned reads (M)", te_m), ("gene-assigned reads (M)", gene_m), ("total counted (M)", tot / 1e6)]
    if qc_dir is not None:
        def star(run, key):
            for line in (qc_dir / f"{run}.Log.final.out").read_text().splitlines():
                if key in line:
                    return float(line.split("|")[1].strip()) / 1e6
        series = [("input reads (M)", pd.Series({r: star(r, "Number of input reads") for r in M.columns})),
                  ("uniquely mapped reads (M)", pd.Series({r: star(r, "Uniquely mapped reads number") for r in M.columns}))] + series
    trows = []
    for lab, x in series:
        t = stats.ttest_ind(x[pert], x[ctrl], equal_var=False)
        ci = t.confidence_interval(0.95)
        trows.append({"study_id": study, "quantity": lab, "control_mean": x[ctrl].mean(),
                      "perturbed_mean": x[pert].mean(), "diff": x[pert].mean() - x[ctrl].mean(),
                      "diff_ci95_low": ci.low, "diff_ci95_high": ci.high, "p_welch": t.pvalue})
        print(f"  {lab}: control {x[ctrl].mean():.2f}, perturbed {x[pert].mean():.2f}, Welch p = {t.pvalue:.2g}")
    if qc_dir is not None:
        pd.DataFrame(trows).to_csv(out.with_name(out.stem + "_read_totals.tsv"), sep="\t", index=False)
    top = F.head(10)
    print(f"  top 5 families carry {top.head(5).frac_of_total_delta.sum():.1%}, top 10 {top.frac_of_total_delta.sum():.1%} of the change")
    print(top.round(4).to_string(index=False))
    print(f"  families with share up: {(F.delta_share_pct_points > 0).sum()} of {len(F)}")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3], Path(sys.argv[4]),
         Path(sys.argv[5]) if len(sys.argv) > 5 else None)
