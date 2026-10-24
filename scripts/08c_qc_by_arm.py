#!/usr/bin/env python3
"""Library QC by arm for every zebrafish contrast (review item: symmetric QC).

Per run: STAR unique-mapping %, multi-mapping %, "too short" %, and fastp
duplication rate. Per study: arm means, perturbed - control difference with
Welch 95% CI and p, Mann-Whitney p, point-biserial r with genotype, BH-FDR
across studies. A study is flagged when its unique-mapping gap matches the
GSE234993 pattern: |difference| >= 5 points and |r| >= 0.7.

Usage: 08c_qc_by_arm.py QC_DIR SAMPLES_TSV OUT_DIR
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def star(path: Path) -> dict:
    d = {}
    for line in path.read_text().splitlines():
        if "|" in line:
            k, v = (x.strip() for x in line.split("|", 1))
            d[k] = v
    pct = lambda k: float(d[k].rstrip("%"))  # noqa: E731
    return {"unique_pct": pct("Uniquely mapped reads %"),
            "multi_pct": pct("% of reads mapped to multiple loci"),
            "too_short_pct": pct("% of reads unmapped: too short")}


def main(qc: Path, samples: Path, out: Path) -> None:
    S = pd.read_csv(samples, sep="\t")
    rows = []
    for r in S.itertuples():
        d = {"run": r.sample_id, "study_id": r.parental_phenotype, "role": r.treatment,
             **star(qc / f"{r.sample_id}.Log.final.out")}
        fj = qc / f"{r.sample_id}.fastp.json"
        d["dup_pct"] = 100 * json.loads(fj.read_text())["duplication"]["rate"] if fj.exists() else np.nan
        rows.append(d)
    R = pd.DataFrame(rows)
    out.mkdir(parents=True, exist_ok=True)
    R.to_csv(out / "08c_qc_per_run.tsv", sep="\t", index=False)

    res = []
    for sid, d in R.groupby("study_id"):
        p, c = d[d.role == "perturbed"], d[d.role == "control"]
        row = {"study_id": sid, "n_perturbed": len(p), "n_control": len(c)}
        for m in ("unique_pct", "dup_pct", "too_short_pct"):
            a, b = p[m].dropna(), c[m].dropna()
            row[f"{m}_perturbed"], row[f"{m}_control"] = a.mean(), b.mean()
            row[f"{m}_diff"] = a.mean() - b.mean()
            if len(a) > 1 and len(b) > 1 and (a.var() > 0 or b.var() > 0):
                w = stats.ttest_ind(a, b, equal_var=False)
                ci = w.confidence_interval(0.95)
                row[f"{m}_diff_ci95_low"], row[f"{m}_diff_ci95_high"], row[f"{m}_p_welch"] = ci.low, ci.high, w.pvalue
                row[f"{m}_p_mannwhitney"] = stats.mannwhitneyu(a, b).pvalue
            if m == "unique_pct":
                g = (d.role == "perturbed").astype(int)
                row["r_genotype_unique"] = stats.pointbiserialr(g, d[m]).statistic
        res.append(row)
    Q = pd.DataFrame(res)
    for m in ("unique_pct", "dup_pct"):
        pv = Q[f"{m}_p_welch"]
        Q[f"{m}_fdr_welch"] = np.nan
        Q.loc[pv.notna(), f"{m}_fdr_welch"] = multipletests(pv.dropna(), method="fdr_bh")[1]
    Q["flag_mapping_gap"] = (Q["unique_pct_diff"].abs() >= 5) & (Q["r_genotype_unique"].abs() >= 0.7)
    Q.to_csv(out / "08c_qc_by_arm.tsv", sep="\t", index=False)
    cols = ["study_id", "n_perturbed", "n_control", "unique_pct_perturbed", "unique_pct_control", "unique_pct_diff",
            "unique_pct_diff_ci95_low", "unique_pct_diff_ci95_high", "unique_pct_p_welch", "unique_pct_fdr_welch",
            "r_genotype_unique", "dup_pct_perturbed", "dup_pct_control", "dup_pct_p_welch", "flag_mapping_gap"]
    print(Q[cols].round(3).to_string(index=False))


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]))
