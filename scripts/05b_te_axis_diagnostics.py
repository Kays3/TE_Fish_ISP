#!/usr/bin/env python3
"""05b — Diagnose the global TE-load axis behind a bulk linear ISP run.

Within the strata used by 05 (default tissue x treatment), takes TE PC1 and
reports, with statistics:
  * variance explained and the share of families loading in the same direction
    (binomial test vs 0.5);
  * correlation of PC1 with gene PC1, TE read fraction and library QC
    (Pearson r, Fisher-z 95% CI, p; degrees of freedom reduced by the number
    of strata because every feature was centred within stratum; QC values
    are standardized within the same strata);
  * correlation of every regulator with PC1, compared with 2,000 random
    genes (Mann-Whitney U on |r|, two-sided) and each regulator's empirical
    percentile among random genes.

Outputs (tables_dir): 05b_axis_summary.tsv, 05b_axis_correlates.tsv,
05b_regulator_axis_corr.tsv

Usage:
  python scripts/05b_te_axis_diagnostics.py --config <run>.yaml [--qc qc.tsv]
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from utils.common import NON_TE_CLASSES, load_config, provenance, setup_logging, write_json  # noqa: E402

_spec = importlib.util.spec_from_file_location("bulk05", HERE / "05_bulk_linear_isp.py")
bulk05 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bulk05)


def corr_stats(a: np.ndarray, b: np.ndarray, n_strata: int) -> dict:
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    r = float(np.corrcoef(a, b)[0, 1])
    n_eff = len(a) - n_strata          # one mean removed per stratum
    z, se = np.arctanh(np.clip(r, -0.999999, 0.999999)), 1 / np.sqrt(n_eff - 3)
    t = r * np.sqrt((n_eff - 2) / max(1e-12, 1 - r ** 2))
    return {"r": r, "ci95_low": float(np.tanh(z - 1.96 * se)), "ci95_high": float(np.tanh(z + 1.96 * se)),
            "p": float(2 * stats.t.sf(abs(t), n_eff - 2)), "n": int(len(a)), "df": int(n_eff - 2)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--config", default=None)
    ap.add_argument("--qc", default=None, help="TSV: sample_id + numeric QC columns")
    ap.add_argument("--n-random", type=int, default=2000)
    args = ap.parse_args()
    cfg = load_config(args.config)
    log = setup_logging(cfg, "05b_te_axis_diagnostics")
    bc = dict(cfg["bulk_isp"])
    bc["counts"] = cfg.path("bulk_isp", "counts")
    tables = cfg.path("paths", "tables_dir")
    tables.mkdir(parents=True, exist_ok=True)

    counts = bulk05.load_counts(bc, log)
    samples = pd.read_csv(cfg.path("bulk_isp", "samples"), sep="\t", dtype=str).set_index(bc["sample_id_col"])
    counts = counts[samples.index.intersection(counts.columns)]
    samples = samples.loc[counts.columns]
    genes, tes = bulk05.split_features(counts, NON_TE_CLASSES | set(cfg.te.get("exclude_classes", [])))
    x = bulk05.normalize(counts, samples, genes, tes, bc, log)
    regs = bulk05.pick_regulators(x.columns, cfg, bc, log)
    strata = bulk05.strata_of(samples.loc[x.index], bc)
    n_strata = int(strata.nunique())
    te_cols = [c for c in x.columns if c in tes]
    gene_cols = [c for c in x.columns if c in genes]

    T = x[te_cols].to_numpy()
    u, sv, vt = np.linalg.svd(T - T.mean(0), full_matrices=False)
    pc1 = u[:, 0] * sv[0]
    if np.median(vt[0]) < 0:
        pc1, vt = -pc1, -vt
    n_pos = int((vt[0] > 0).sum())
    G = x[gene_cols].to_numpy()
    ug, sg, _ = np.linalg.svd(G - G.mean(0), full_matrices=False)
    summary = pd.DataFrame([{
        "n_samples": len(x), "n_strata": n_strata, "n_te_families": len(te_cols),
        "te_pc1_var_explained": float(sv[0] ** 2 / (sv ** 2).sum()),
        "te_pc2_var_explained": float(sv[1] ** 2 / (sv ** 2).sum()),
        "frac_families_same_sign": n_pos / len(te_cols),
        "p_binom_same_sign": float(stats.binomtest(n_pos, len(te_cols), 0.5).pvalue),
        "gene_pc1_var_explained": float(sg[0] ** 2 / (sg ** 2).sum()),
    }])
    summary.to_csv(tables / "05b_axis_summary.tsv", sep="\t", index=False)

    def within(v: pd.Series) -> np.ndarray:
        v = pd.to_numeric(v, errors="coerce")
        return v.groupby(strata).transform(lambda a: (a - a.mean()) / (a.std() or 1)).to_numpy()

    rows = [{"variable": "gene PC1", **corr_stats(pc1, ug[:, 0], n_strata)}]
    rows[0]["r"] = abs(rows[0]["r"])          # PC sign is arbitrary
    rows[0]["ci95_low"], rows[0]["ci95_high"] = sorted(map(abs, (rows[0]["ci95_low"], rows[0]["ci95_high"])))
    rows.append({"variable": "TE read fraction",
                 **corr_stats(pc1, within(counts.loc[tes].sum() / counts.sum()), n_strata)})
    rows.append({"variable": "library size (counted reads)", **corr_stats(pc1, within(counts.sum()), n_strata)})
    if args.qc:
        qc = pd.read_csv(args.qc, sep="\t").set_index("sample_id").reindex(x.index)
        for c in qc.columns:
            if pd.api.types.is_numeric_dtype(qc[c]):
                rows.append({"variable": c, **corr_stats(pc1, within(qc[c]), n_strata)})
    cor = pd.DataFrame(rows)
    cor.to_csv(tables / "05b_axis_correlates.tsv", sep="\t", index=False)

    rng = np.random.default_rng(cfg.project["seed"])
    pool = [g for g in gene_cols if g not in regs]
    rand = rng.choice(pool, min(args.n_random, len(pool)), replace=False)
    r_rand = np.array([np.corrcoef(pc1, x[g])[0, 1] for g in rand])
    reg = pd.DataFrame([{"regulator": r, **corr_stats(pc1, x[r].to_numpy(), n_strata)} for r in regs])
    reg["pct_random_abs_r_ge"] = [float((np.abs(r_rand) >= abs(v)).mean()) for v in reg["r"]]
    reg = reg.sort_values("r")
    reg.to_csv(tables / "05b_regulator_axis_corr.tsv", sep="\t", index=False)
    mw = stats.mannwhitneyu(np.abs(reg["r"]), np.abs(r_rand), alternative="two-sided")
    write_json(tables / "05b_provenance.json", provenance(
        cfg, regulators_vs_random_abs_r={
            "median_abs_r_regulators": float(np.median(np.abs(reg["r"]))),
            "median_abs_r_random": float(np.median(np.abs(r_rand))),
            "mannwhitney_U": float(mw.statistic), "p": float(mw.pvalue),
            "n_regulators": len(reg), "n_random": len(r_rand)}))
    log.info("TE PC1 %.1f%%; |r| regulators median %.2f vs random %.2f (Mann-Whitney p=%.3g)",
             100 * summary.te_pc1_var_explained.iat[0], np.median(np.abs(reg["r"])),
             np.median(np.abs(r_rand)), mw.pvalue)


if __name__ == "__main__":
    main()
