#!/usr/bin/env python3
"""07 — Which genes respond to treatment, and do they account for the
treatment's effect on global TE load?

For each tissue, each candidate gene set (config/gene_sets_*.yaml) and each
treatment contrast (developmental / intergenerational vs control):

  A. Treatment effect on every gene: log2 expression ~ family + treatment
     (frozen OA design, family = parental line nested in phenotype), giving
     log2FC, 95% CI, p and BH-FDR across the candidate genes. Set-level
     shift is tested competitively: the mean t of the set against 10,000
     random gene sets matched gene-for-gene on expression level.
  B. Within-group association of each gene with the TE axis:
     axis ~ family + treatment + gene, giving the partial slope (per SD of
     the gene), p and FDR.
  C. Mediation (difference in coefficients). Total treatment effect on the
     axis c; with the gene added, c' = c - beta*a (beta = B slope, a =
     A effect, exact OLS identity). The indirect share (beta*a)/c gets a
     bootstrap 95% CI (fish resampled within family x treatment cells) and
     an empirical p against expression-matched random genes. Set-level
     mediation adds the set's top principal components together and is
     compared with random sets of the same size.

Because the TE axis correlates with most of the transcriptome (05b), every
gene-level and set-level claim is made relative to random genes, never
against zero alone. Mediation here is statistical, not causal: it assumes no
unmeasured confounding of gene and TE axis, which the global-axis structure
makes doubtful.

Usage:
  python scripts/07_treatment_gene_axis.py --config runs/<run>.yaml \
      --gene-sets config/gene_sets_apoly.yaml [--n-random-sets 10000] [--n-boot 1000]
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from utils.common import NON_TE_CLASSES, load_config, provenance, setup_logging, write_json  # noqa: E402


def _load(name: str, file: str):
    spec = importlib.util.spec_from_file_location(name, HERE / file)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


bulk05 = _load("bulk05", "05_bulk_linear_isp.py")
axis06 = _load("axis06", "06_te_axis_treatment.py")

CONTRASTS = {"developmental": "developmental_vs_control", "intergenerational": "intergenerational_vs_control"}
SET_COLORS = {"te_regulators": "#4d4d4d", "gaba": "#2166ac", "thyroid": "#b2182b", "glut": "#1b7837"}


def bh(p):
    p = np.asarray(p, float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if ok.any():
        out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out


def design(meta: pd.DataFrame) -> tuple[np.ndarray, dict]:
    fam = pd.get_dummies(meta["fam_id"], drop_first=True, dtype=float)
    trt = pd.get_dummies(pd.Categorical(meta["treatment"], ["control", *CONTRASTS]), dtype=float)
    X = np.column_stack([np.ones(len(meta)), fam.to_numpy(), trt[list(CONTRASTS)].to_numpy()])
    idx = {lv: X.shape[1] - len(CONTRASTS) + i for i, lv in enumerate(CONTRASTS)}
    return X, idx


def ols_many(X: np.ndarray, Y: np.ndarray):
    """Coefficients, SEs and residual df for many outcomes at once."""
    XtX_inv = np.linalg.pinv(X.T @ X)
    B = XtX_inv @ X.T @ Y
    R = Y - X @ B
    df = X.shape[0] - np.linalg.matrix_rank(X)
    s2 = (R ** 2).sum(0) / df
    se = np.sqrt(np.outer(np.diag(XtX_inv), s2))
    return B, se, df


def partial_slopes(X: np.ndarray, y: np.ndarray, G: np.ndarray):
    """Slope of y on each column of G given X (Frisch-Waugh), with t-test."""
    M = np.eye(len(y)) - X @ np.linalg.pinv(X)
    ry, RG = M @ y, M @ G
    ss = (RG ** 2).sum(0)
    beta = RG.T @ ry / np.where(ss > 0, ss, np.nan)
    df = len(y) - np.linalg.matrix_rank(X) - 1
    resid = ry[:, None] - RG * beta
    se = np.sqrt((resid ** 2).sum(0) / df / ss)
    return beta, se, df


def matched_draw(level: pd.Series, targets: list[str], pool: list[str], rng, n_bins=20) -> np.ndarray:
    """For each target gene, draw a pool gene from the same expression bin."""
    bins = pd.qcut(level.rank(method="first"), n_bins, labels=False)
    by_bin = {b: np.array([g for g in pool if bins[g] == b]) for b in range(n_bins)}
    return np.array([rng.choice(by_bin[bins[t]]) for t in targets])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--config", default=None)
    ap.add_argument("--gene-sets", required=True)
    ap.add_argument("--n-random-sets", type=int, default=10000)
    ap.add_argument("--n-random-genes", type=int, default=3000)
    ap.add_argument("--n-boot", type=int, default=1000)
    args = ap.parse_args()
    cfg = load_config(args.config)
    log = setup_logging(cfg, "07_treatment_gene_axis")
    bc = dict(cfg["bulk_isp"])
    bc["counts"] = cfg.path("bulk_isp", "counts")
    tables, figs = cfg.path("paths", "tables_dir"), cfg.path("paths", "figures_dir")
    tables.mkdir(parents=True, exist_ok=True)
    figs.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(cfg.project["seed"])
    gsets = yaml.safe_load(open(args.gene_sets))

    counts = bulk05.load_counts(bc, log)
    samples = pd.read_csv(cfg.path("bulk_isp", "samples"), sep="\t", dtype=str).set_index(bc["sample_id_col"])
    counts = counts[samples.index.intersection(counts.columns)]
    samples = samples.loc[counts.columns]
    genes, tes = bulk05.split_features(counts, NON_TE_CLASSES | set(cfg.te.get("exclude_classes", [])))
    sc, _, _ = axis06.axis_scores(counts, samples, genes, tes, bc, log)
    expr = bulk05.normalize(counts, samples, genes, tes, dict(bc, stratify_by=["tissue"], scale=False), log)
    gene_cols = [c for c in expr.columns if c in genes]
    meta = samples.assign(fam_id=samples["parental_phenotype"] + "_" + samples["family"].astype(str))

    member = {g: s for s, spec in gsets.items() for g in spec["genes"]}
    missing = sorted(g for g in member if g not in gene_cols)
    if missing:
        log.warning("Set genes absent after global filtering: %s", ", ".join(missing))

    eff_rows, assoc_rows, med_rows, setshift_rows, setmed_rows = [], [], [], [], []
    for tissue in sorted(meta["tissue"].unique()):
        idx_t = meta.index[meta["tissue"] == tissue]
        m = meta.loc[idx_t]
        raw_mean = counts.loc[gene_cols, idx_t].mean(axis=1)
        expressed = [g for g in gene_cols if raw_mean[g] >= bc["min_count"]]
        E = expr.loc[idx_t, expressed]
        Ez = (E - E.mean()) / E.std(ddof=1).replace(0, np.nan)
        Ez = Ez.dropna(axis=1)
        expressed = list(Ez.columns)
        level = np.log2(raw_mean[expressed] + 1)
        X, tix = design(m)
        y = sc.loc[idx_t, "axis_score"].to_numpy()
        set_genes = [g for g in expressed if g in member]
        pool = [g for g in expressed if g not in member]
        pos = {g: i for i, g in enumerate(expressed)}
        bins = pd.qcut(level.rank(method="first"), 20, labels=False)
        pool_bins = bins[pool]
        pool_by_bin = {b: pool_bins.index[pool_bins == b].to_numpy() for b in range(20)}
        log.info("%s: %d expressed genes; %d candidate genes in sets", tissue, len(expressed), len(set_genes))

        # A. treatment effects on all expressed genes (log2 scale, unstandardized)
        B, SE, df = ols_many(X, E[expressed].to_numpy())
        tcrit = stats.t.ppf(0.975, df)
        T = {lv: B[j] / SE[j] for lv, j in tix.items()}
        for lv, j in tix.items():
            pv = 2 * stats.t.sf(np.abs(T[lv]), df)
            for k, g in enumerate(expressed):
                if g in member:
                    eff_rows.append({"tissue": tissue, "contrast": CONTRASTS[lv], "gene": g, "set": member[g],
                                     "log2fc": B[j, k], "ci95_low": B[j, k] - tcrit * SE[j, k],
                                     "ci95_high": B[j, k] + tcrit * SE[j, k], "t": T[lv][k], "p": pv[k],
                                     "mean_count": float(raw_mean[g])})
            # competitive set test: mean t vs expression-matched random sets
            tser = pd.Series(T[lv], index=expressed)
            for sname, spec in gsets.items():
                sg = [g for g in spec["genes"] if g in tser.index]
                if len(sg) < 3:
                    continue
                obs = float(tser[sg].mean())
                obs_abs = float(tser[sg].abs().mean())
                by_bin = {b: tser[pool_by_bin[b]].to_numpy() for b in range(20)}
                draws = np.column_stack([rng.choice(by_bin[bins[g]], args.n_random_sets) for g in sg])
                null, null_abs = draws.mean(1), np.abs(draws).mean(1)
                setshift_rows.append({
                    "tissue": tissue, "contrast": CONTRASTS[lv], "set": sname, "n_genes": len(sg),
                    "mean_t": obs, "null_mean": float(null.mean()), "null_sd": float(null.std()),
                    "z_vs_null": (obs - null.mean()) / null.std(),
                    "p_directional": float((1 + (np.abs(null - null.mean()) >= abs(obs - null.mean())).sum())
                                           / (1 + len(null))),
                    "mean_abs_t": obs_abs, "null_mean_abs_t": float(null_abs.mean()),
                    "p_any_direction": float((1 + (null_abs >= obs_abs).sum()) / (1 + len(null_abs))),
                    "frac_up": float((tser[sg] > 0).mean()),
                    "p_binom_up": float(stats.binomtest(int((tser[sg] > 0).sum()), len(sg)).pvalue)})

        # B. within-group partial association with the TE axis (per SD of gene)
        Gz = Ez[expressed].to_numpy()
        beta, bse, bdf = partial_slopes(X, y, Gz)
        bp = 2 * stats.t.sf(np.abs(beta / bse), bdf)
        bser = pd.Series(beta, index=expressed)
        for k, g in enumerate(expressed):
            if g in member:
                q = stats.t.ppf(0.975, bdf)
                assoc_rows.append({"tissue": tissue, "gene": g, "set": member[g], "beta_per_sd": beta[k],
                                   "ci95_low": beta[k] - q * bse[k], "ci95_high": beta[k] + q * bse[k],
                                   "p": bp[k], "pct_random_abs_ge": np.nan})
        rand_genes = rng.choice(pool, min(args.n_random_genes, len(pool)), replace=False)
        rabs = np.abs(bser[rand_genes].to_numpy())
        for r in assoc_rows:
            if r["tissue"] == tissue:
                k_ge = int((rabs >= abs(r["beta_per_sd"])).sum())
                r["pct_random_abs_ge"] = k_ge / len(rabs)
                r["p_emp_vs_random"] = (1 + k_ge) / (1 + len(rabs))

        # C. mediation of each contrast's effect on the axis
        Bz, _, _ = ols_many(X, Gz)                       # a: treatment effect on standardized gene
        c_tot, c_se, c_df = ols_many(X, y[:, None])
        cells = (m["fam_id"] + "|" + m["treatment"]).to_numpy()
        cell_idx = [np.flatnonzero(cells == c) for c in np.unique(cells)]
        sg_idx = [pos[g] for g in set_genes]
        boot_ind = {lv: np.empty((args.n_boot, len(sg_idx))) for lv in tix}
        boot_c = {lv: np.empty(args.n_boot) for lv in tix}
        for b in range(args.n_boot):
            ii = np.concatenate([rng.choice(ci, len(ci)) for ci in cell_idx])
            Xb, yb, Gb = X[ii], y[ii], Gz[ii][:, sg_idx]
            cb, _, _ = ols_many(Xb, yb[:, None])
            ab, _, _ = ols_many(Xb, Gb)
            betab, _, _ = partial_slopes(Xb, yb, Gb)
            for lv, j in tix.items():
                boot_ind[lv][b] = betab * ab[j]
                boot_c[lv][b] = cb[j, 0]
        for lv, j in tix.items():
            c = float(c_tot[j, 0])
            c_p = float(2 * stats.t.sf(abs(c / c_se[j, 0]), c_df))
            ind_all = pd.Series(beta * Bz[j], index=expressed)
            for s_i, g in enumerate(set_genes):
                ind = float(ind_all[g])
                same_bin = [h for h in rand_genes if bins[h] == bins[g]]
                nullv = ind_all[same_bin].to_numpy()
                bi = boot_ind[lv][:, s_i]
                prop_boot = bi / boot_c[lv]
                med_rows.append({
                    "tissue": tissue, "contrast": CONTRASTS[lv], "gene": g, "set": member[g],
                    "total_effect": c, "total_effect_p": c_p, "a_trt_on_gene_sd": float(Bz[j, pos[g]]),
                    "b_gene_on_axis": float(bser[g]), "indirect": ind,
                    "indirect_ci95_low": float(np.quantile(bi, 0.025)),
                    "indirect_ci95_high": float(np.quantile(bi, 0.975)),
                    "prop_mediated": ind / c if c else np.nan,
                    "prop_ci95_low": float(np.quantile(prop_boot, 0.025)),
                    "prop_ci95_high": float(np.quantile(prop_boot, 0.975)),
                    "n_null_same_bin": len(nullv),
                    "null_indirect_q025": float(np.quantile(nullv, 0.025)),
                    "null_indirect_q975": float(np.quantile(nullv, 0.975)),
                    "p_emp_vs_random": float((1 + (np.abs(nullv) >= abs(ind)).sum()) / (1 + len(nullv)))})

            # set-level mediation: add the set's top-k PCs together
            def attenuation(cols: list[int]) -> float:
                G = Gz[:, cols]
                u, s_, _ = np.linalg.svd(G - G.mean(0), full_matrices=False)
                k = min(3, len(cols))
                Xm = np.column_stack([X, u[:, :k]])
                cm, _, _ = ols_many(Xm, y[:, None])
                return (c - float(cm[j, 0])) / c if c else np.nan
            for sname, spec in gsets.items():
                sg = [pos[g] for g in spec["genes"] if g in pos]
                if len(sg) < 3:
                    continue
                obs = attenuation(sg)
                nulls = np.empty(1000)
                for r_ in range(1000):
                    draw = [rng.choice(pool_by_bin[bins[expressed[i]]]) for i in sg]
                    nulls[r_] = attenuation([pos[d] for d in draw])
                setmed_rows.append({"tissue": tissue, "contrast": CONTRASTS[lv], "set": sname, "n_genes": len(sg),
                                    "total_effect": c, "total_effect_p": c_p, "prop_mediated_set_pcs": obs,
                                    "null_median": float(np.median(nulls)),
                                    "null_q975": float(np.quantile(nulls, 0.975)),
                                    "p_emp_vs_random_sets": float((1 + (nulls >= obs).sum()) / (1 + len(nulls)))})
        log.info("%s done", tissue)

    eff = pd.DataFrame(eff_rows)
    eff["fdr"] = np.nan
    for _, ix in eff.groupby(["tissue", "contrast"]).groups.items():
        eff.loc[ix, "fdr"] = bh(eff.loc[ix, "p"])
    assoc = pd.DataFrame(assoc_rows)
    assoc["fdr"] = np.nan
    assoc["fdr_emp"] = np.nan
    for _, ix in assoc.groupby("tissue").groups.items():
        assoc.loc[ix, "fdr"] = bh(assoc.loc[ix, "p"])
        assoc.loc[ix, "fdr_emp"] = bh(assoc.loc[ix, "p_emp_vs_random"])
    med = pd.DataFrame(med_rows)
    med["fdr_emp"] = np.nan
    for _, ix in med.groupby(["tissue", "contrast"]).groups.items():
        med.loc[ix, "fdr_emp"] = bh(med.loc[ix, "p_emp_vs_random"])
    shift = pd.DataFrame(setshift_rows)
    for col in ("p_directional", "p_any_direction"):
        shift[col.replace("p_", "fdr_")] = bh(shift[col])
    setmed = pd.DataFrame(setmed_rows)
    # A proportion mediated is only interpretable when there is a total effect
    # to mediate; FDR is computed over those rows only.
    setmed["interpretable"] = setmed["total_effect_p"] < 0.05
    setmed["fdr"] = np.nan
    setmed.loc[setmed.interpretable, "fdr"] = bh(setmed.loc[setmed.interpretable, "p_emp_vs_random_sets"])
    med["interpretable"] = med["total_effect_p"] < 0.05
    med.loc[~med.interpretable, "fdr_emp"] = np.nan
    for _, ix in med[med.interpretable].groupby(["tissue", "contrast"]).groups.items():
        med.loc[ix, "fdr_emp"] = bh(med.loc[ix, "p_emp_vs_random"])
    for name, df_ in (("07_gene_treatment_effects", eff), ("07_set_treatment_shift", shift),
                      ("07_gene_axis_association", assoc), ("07_gene_mediation", med),
                      ("07_set_mediation", setmed)):
        df_.to_csv(tables / f"{name}.tsv", sep="\t", index=False)
    make_figures(eff, shift, assoc, med, setmed, gsets, figs)
    write_json(tables / "07_provenance.json", provenance(
        cfg, gene_sets=args.gene_sets, n_random_sets=args.n_random_sets, n_boot=args.n_boot,
        n_random_genes=args.n_random_genes, missing_genes=missing))
    log.info("Wrote 07 tables and figures")


# ---------------------------------------------------------------------------
def make_figures(eff, shift, assoc, med, setmed, gsets, figs):  # noqa: C901
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                         "font.family": "DejaVu Sans"})
    cols = [(t, c) for t in ("brain", "gill", "liver") for c in CONTRASTS.values()]
    col_lab = [f"{t}\n{'dev' if c.startswith('dev') else 'intergen'}" for t, c in cols]

    # Figure 1: heatmap of treatment effects (log2FC) with significance marks
    order = [(s, g) for s, spec in gsets.items() for g in spec["genes"]]
    order = [(s, g) for s, g in order if g in set(eff["gene"])]
    M = np.full((len(order), len(cols)), np.nan)
    S = np.full(M.shape, "", dtype=object)
    for i, (s, g) in enumerate(order):
        for j, (t, c) in enumerate(cols):
            r = eff[(eff.gene == g) & (eff.tissue == t) & (eff.contrast == c)]
            if len(r):
                M[i, j] = r.log2fc.iat[0]
                S[i, j] = "**" if r.fdr.iat[0] < 0.05 else ("·" if r.p.iat[0] < 0.05 else "")
    lim = np.nanquantile(np.abs(M), 0.98)
    fig, ax = plt.subplots(figsize=(4.8, 0.13 * len(order) + 1.6))
    im = ax.imshow(M, aspect="auto", cmap="RdBu_r", norm=TwoSlopeNorm(0, -lim, lim), interpolation="none")
    for i in range(len(order)):
        for j in range(len(cols)):
            if S[i, j]:
                ax.text(j, i, S[i, j], ha="center", va="center", fontsize=6)
    ax.set_xticks(range(len(cols)), col_lab, fontsize=7)
    ax.set_yticks(range(len(order)), [g for _, g in order], fontsize=5.5)
    for tick, (s, _) in zip(ax.get_yticklabels(), order):
        tick.set_color(SET_COLORS.get(s, "k"))
    bounds = np.cumsum([sum(1 for s2, _ in order if s2 == s) for s in gsets])[:-1]
    for b in bounds:
        ax.axhline(b - 0.5, color="k", lw=0.6)
    for b in (1.5, 3.5):
        ax.axvline(b, color="k", lw=0.6)
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("log2 fold change vs control")
    ax.set_title("Treatment effects on candidate genes\n** BH-FDR < 0.05 (per tissue × contrast), · p < 0.05",
                 fontsize=8)
    fig.tight_layout()
    fig.savefig(figs / "07_fig1_gene_effects_heatmap.png", dpi=300)
    fig.savefig(figs / "07_fig1_gene_effects_heatmap.pdf")
    plt.close(fig)

    # Figure 2: set-level competitive shift (z vs matched random sets)
    fig, axes = plt.subplots(1, len(gsets), figsize=(9.5, 2.9), sharey=True)
    ylab = [f"{t} {'dev' if c.startswith('dev') else 'intergen'}" for t, c in cols]
    for ax, s in zip(axes, gsets):
        d = shift[shift.set == s].set_index(["tissue", "contrast"])
        yy = np.arange(len(cols))
        ax.axvspan(-1.96, 1.96, color="0.92", zorder=0)
        for k, (t, c) in enumerate(cols):
            if (t, c) not in d.index:
                continue
            r = d.loc[(t, c)]
            ax.barh(k, r.z_vs_null, color=SET_COLORS[s], height=0.65)
            ax.text(1.02, k, f"p={r.p_directional:.2g}" + ("*" if r.fdr_directional < 0.05 else ""),
                    transform=ax.get_yaxis_transform(), va="center", fontsize=6, clip_on=False)
        ax.axvline(0, color="k", lw=0.6)
        n = int(d.n_genes.max()) if len(d) else 0
        ax.set_title(f"{s} (n≤{n})", fontsize=8)
        ax.set_xlabel("z vs matched random sets", fontsize=7)
        ax.set_yticks(yy, ylab, fontsize=7)
        ax.set_ylim(len(cols) - 0.5, -0.5)
    fig.suptitle("Gene-set shift with treatment: set mean t vs 10,000 expression-matched random sets\n"
                 "grey band = |z| < 1.96; raw p shown, * BH-FDR < 0.05 across all set × tissue × contrast tests",
                 fontsize=8)
    fig.tight_layout(rect=(0, 0, 0.97, 1))
    fig.subplots_adjust(wspace=0.55)
    fig.savefig(figs / "07_fig2_set_shift.png", dpi=300)
    fig.savefig(figs / "07_fig2_set_shift.pdf")
    plt.close(fig)

    # Figure 3: gene-level mediation for contrasts with a nominal total effect
    targets = med.drop_duplicates(["tissue", "contrast"])
    targets = targets[targets.total_effect_p < 0.05][["tissue", "contrast"]].values.tolist()
    if targets:
        fig, axes = plt.subplots(1, len(targets), figsize=(max(5.8, 3.8 * len(targets)), 0.11 * len(order) + 2.0),
                                 squeeze=False)
        for ax, (t, c) in zip(axes[0], targets):
            d = med[(med.tissue == t) & (med.contrast == c)].copy()
            d["key"] = pd.Categorical(d.gene, [g for _, g in order], ordered=True)
            d = d.sort_values("key")
            yy = np.arange(len(d))
            ax.hlines(yy, d.null_indirect_q025, d.null_indirect_q975, color="0.82", lw=4, zorder=0,
                      label="95% range, matched random genes")
            ax.errorbar(d.indirect, yy, xerr=[d.indirect - d.indirect_ci95_low, d.indirect_ci95_high - d.indirect],
                        fmt="none", ecolor="0.4", lw=0.7, capsize=1.5)
            ax.scatter(d.indirect, yy, c=[SET_COLORS[s] for s in d.set], s=10, zorder=3)
            for k, r in enumerate(d.itertuples()):
                if r.fdr_emp < 0.05:
                    ax.text(r.indirect_ci95_high, k, " *", va="center", fontsize=7)
            ax.axvline(0, color="k", lw=0.6)
            ax.set_yticks(yy, d.gene, fontsize=5.5)
            ax.invert_yaxis()
            tot = d.total_effect.iat[0]
            ax.set_xlabel("indirect effect on TE axis (SD) = a × b")
            ax.set_title(f"{t}, {c.split('_')[0]} (total effect {tot:+.2f} SD)", fontsize=8)
        axes[0][0].legend(frameon=False, fontsize=6, loc="lower right")
        fig.suptitle("Gene-level mediation of the treatment effect on the global TE axis\n"
                     "point = estimate; bar = bootstrap 95% CI; grey = 95% range of matched random genes\n"
                     "* BH-FDR < 0.05 vs matched random genes", fontsize=7.5)
        fig.tight_layout()
        fig.savefig(figs / "07_fig3_gene_mediation.png", dpi=300)
        fig.savefig(figs / "07_fig3_gene_mediation.pdf")
        plt.close(fig)

    # Figure 5: set-level mediation where a total effect exists
    sm = setmed[setmed.interpretable]
    if len(sm):
        fig, ax = plt.subplots(figsize=(3.6, 2.4))
        for k, r in enumerate(sm.itertuples()):
            ax.plot([r.null_median, r.null_q975], [k, k], color="0.75", lw=5, solid_capstyle="butt", zorder=0)
            ax.scatter(r.prop_mediated_set_pcs, k, color=SET_COLORS[r.set], s=30, zorder=3)
            ax.text(max(r.prop_mediated_set_pcs, r.null_q975) + 0.05, k,
                    f"p={r.p_emp_vs_random_sets:.2g}, FDR={r.fdr:.2g}", va="center", fontsize=6)
        ax.set_yticks(range(len(sm)), [f"{r.tissue} {r.contrast.split('_')[0][:3]}: {r.set}" for r in sm.itertuples()],
                      fontsize=6.5)
        ax.axvline(0, color="k", lw=0.6)
        ax.set_xlabel("share of treatment effect absorbed by set (top PCs)")
        ax.set_title("Set-level mediation vs matched random sets\n(grey: random-set median to 97.5th percentile)",
                     fontsize=7.5)
        ax.invert_yaxis()
        fig.tight_layout()
        fig.savefig(figs / "07_fig5_set_mediation.png", dpi=300)
        fig.savefig(figs / "07_fig5_set_mediation.pdf")
        plt.close(fig)

    # Figure 4: within-group association with the TE axis, by set
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.6), sharey=True)
    for ax, t in zip(axes, ("brain", "gill", "liver")):
        d = assoc[assoc.tissue == t]
        for k, s in enumerate(gsets):
            v = d[d.set == s]
            x = k + np.random.default_rng(k).uniform(-0.18, 0.18, len(v))
            ax.scatter(x, v.beta_per_sd, s=9, c=SET_COLORS[s], alpha=0.85)
            sig = v.fdr_emp < 0.05
            ax.scatter(x[sig.to_numpy()], v.beta_per_sd[sig], s=28, facecolors="none", edgecolors="k", lw=0.6)
        ax.axhline(0, color="k", lw=0.6)
        ax.set_xticks(range(len(gsets)), list(gsets), rotation=30, fontsize=6.5)
        ax.set_title(t, fontsize=8)
    axes[0].set_ylabel("partial slope of TE axis on gene (SD per SD)")
    fig.suptitle("Within-group association with the global TE axis\n"
                 "(ringed: BH-FDR < 0.05 vs 3,000 random genes; none in this run if no rings)", fontsize=8)
    fig.tight_layout()
    fig.savefig(figs / "07_fig4_axis_association.png", dpi=300)
    fig.savefig(figs / "07_fig4_axis_association.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
