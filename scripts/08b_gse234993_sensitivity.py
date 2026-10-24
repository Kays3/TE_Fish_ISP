#!/usr/bin/env python3
"""GSE234993 (ZF_uhrf1_liver_2023) sensitivity table (DC-16, DC-18).

Applies the 08 family rules (strict TE classes, baseMean >= 10; called at
padj < 0.05 and |log2FC| >= 0.5) to six DESeq2 runs of this study: paired-end
and read-1-only alignments, each with ~ block + role (primary), ~ block + uniq +
role and ~ uniq + role (uniq = unique-mapping %, post hoc). Also writes the
per-library mapping table with the genotype tests.

Usage: 08b_gse234993_sensitivity.py IN_DIR OUT_DIR
IN_DIR holds te_de/, r1only/, mapcov_pe/, mapcov_r1only/ as written on the DESeq2 host.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils.common import NON_TE_CLASSES  # noqa: E402

S = "ZF_uhrf1_liver_2023"
RUNS = [  # (alignment, design, path relative to IN_DIR)
    ("paired", "~ block + role", f"te_de/{S}.tsv"),
    ("paired", "~ block + uniq + role", f"mapcov_pe/{S}.pe.mapcov_block.tsv"),
    ("paired", "~ uniq + role", f"mapcov_pe/{S}.pe.mapcov_noblock.tsv"),
    ("read 1 only", "~ block + role", f"r1only/{S}.r1only.tsv"),
    ("read 1 only", "~ block + uniq + role", f"mapcov_r1only/{S}.r1only.mapcov_block.tsv"),
    ("read 1 only", "~ uniq + role", f"mapcov_r1only/{S}.r1only.mapcov_noblock.tsv"),
]


def summarise(path: Path) -> dict:
    de = pd.read_csv(path, sep="\t")
    parts = de["feature"].astype(str).str.split(":")
    cls = parts.str[2].where(parts.str.len() >= 3).str.rstrip("?")
    strict = cls.notna() & ~cls.isin(NON_TE_CLASSES | {"Satellite", "Low_complexity"})
    te = de[strict & (de["baseMean"] >= 10)].dropna(subset=["log2FoldChange"])
    sig = te[(te["padj"] < 0.05) & (te["log2FoldChange"].abs() >= 0.5)]
    up, down = int((sig["log2FoldChange"] > 0).sum()), int((sig["log2FoldChange"] < 0).sum())
    bt = stats.binomtest(up, up + down) if up + down else None
    ci = bt.proportion_ci() if bt else None
    genes = de[~parts.str.len().ge(3)]
    return {"n_te_tested": len(te), "n_sig_up": up, "n_sig_down": down,
            "frac_up": up / (up + down) if up + down else np.nan,
            "frac_up_ci95_low": ci.low if ci else np.nan, "frac_up_ci95_high": ci.high if ci else np.nan,
            "p_binom_up_vs_down": bt.pvalue if bt else np.nan,
            "median_te_lfc": float(te["log2FoldChange"].median()),
            "n_genes_padj05": int((genes["padj"] < 0.05).sum())}


def main(inp: Path, out: Path) -> None:
    rows = [{"study_id": S, "alignment": a, "design": d, "post_hoc": "uniq" in d, **summarise(inp / p)}
            for a, d, p in RUNS]
    T = pd.DataFrame(rows)
    T["fdr_bh_binom"] = multipletests(T["p_binom_up_vs_down"], method="fdr_bh")[1]
    out.mkdir(parents=True, exist_ok=True)
    T.to_csv(out / "08s_gse234993_sensitivity.tsv", sep="\t", index=False)

    m = []
    for a, sub in (("paired", "mapcov_pe"), ("read 1 only", "mapcov_r1only")):
        sh = pd.read_csv(inp / sub / "sheet.tsv", sep="\t").assign(alignment=a)
        m.append(sh)
    M = pd.concat(m).rename(columns={"treatment": "role", "uniq_pct": "unique_mapped_pct"})
    M.to_csv(out / "08s_gse234993_mapping.tsv", sep="\t", index=False)
    tests = []
    for a, d in M.groupby("alignment"):
        p, c = d.loc[d.role == "perturbed", "unique_mapped_pct"], d.loc[d.role == "control", "unique_mapped_pct"]
        w = stats.ttest_ind(p, c, equal_var=False)
        ci = w.confidence_interval(0.95)
        r = stats.pointbiserialr((d.role == "perturbed").astype(int), d.unique_mapped_pct)
        tests.append({"alignment": a, "n_perturbed": len(p), "n_control": len(c),
                      "mean_perturbed": p.mean(), "mean_control": c.mean(), "diff_points": p.mean() - c.mean(),
                      "diff_ci95_low": ci.low, "diff_ci95_high": ci.high, "p_welch": w.pvalue,
                      "p_mannwhitney": stats.mannwhitneyu(p, c).pvalue,
                      "r_genotype_mapping": r.statistic, "p_r": r.pvalue})
    pd.DataFrame(tests).to_csv(out / "08s_gse234993_mapping_tests.tsv", sep="\t", index=False)
    print(T.round(4).to_string(index=False))
    print(pd.DataFrame(tests).round(4).to_string(index=False))


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
