#!/usr/bin/env python3
"""08 — TE response to genetic loss or gain of TE regulators in zebrafish, as
ground truth for the bulk ISP predictions from the A. polyacanthus cohort.

Inputs (produced by pipelines/zebrafish_validation/):
  --te-de     per-study DESeq2 tables (feature, baseMean, log2FoldChange, padj),
              perturbed vs control, apeglm-shrunken (scripts/04_bulk_te_deseq.R)
  --counts    TEcount *.cntTable per run (genes + TEs, --mode multi)
  --samples   samples_all.tsv (sample_id, treatment = control|perturbed,
              parental_phenotype = study_id)
  --manifest  config/zebrafish_perturbation_studies.tsv
  --gene-names gene_id -> gene_name map from the Ensembl GRCz11 GTF (optional;
              enables the perturbed-gene check)
  --fish-05c  directory with 05c_{plain,rmTEpc1}_per_regulator.tsv (optional)

Per study, with statistics:
  * TE families tested (strict TE classes, baseMean >= 10); up / down at
    padj < 0.05 and |log2FC| >= 0.5; exact binomial test of up vs down and the
    95% CI of the fraction up.
  * Median family log2FC with a bootstrap 95% CI over families and a Wilcoxon
    signed-rank p (families are not independent; treat as descriptive).
  * Per TE class: median log2FC, Wilcoxon p, BH across all study x class tests.
  * Global TE share of counted reads per library: logit difference perturbed vs
    control (Welch t, 95% CI), back-transformed to % relative change; BH across
    studies.
  * Perturbed-gene check: the targeted gene's own log2FC (mutants with NMD and
    knockdowns should fall; morpholinos that block translation, and transgenes
    of the human gene, need not change zebrafish mRNA).

Across studies: does the observed direction (TE share up or down, called only
when BH-FDR < 0.05 and |change| >= --min-share-change %) match
(a) the expectation that losing a silencer derepresses TEs, and
(b) the direction predicted for the matching fish regulator by
05_bulk_linear_isp.py (fraction of significant TE families up, plain and
axis-removed runs)? Exact binomial tests on the number of agreements.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils.common import NON_TE_CLASSES  # noqa: E402

FISH_REG = {"DNMT1": ["dnmt1"], "UHRF1": ["uhrf1"], "EZH2": ["ezh2"], "KDM1A": ["kdm1a"], "ATRX": ["atrx"],
            "METTL3": ["mettl3"], "DNMT3B": ["dnmt3bb.1"]}
ZF_GENE = {"DNMT1": ["dnmt1"], "UHRF1": ["uhrf1"], "EZH2": ["ezh2"], "KDM1A": ["kdm1a"], "ATRX": ["atrx"],
           "METTL3": ["mettl3"], "DNMT3B": ["dnmt3bb.1", "dnmt3ba", "dnmt3bb.2", "dnmt3bb.3"]}
CLASS_ORDER = ["LINE", "SINE", "LTR", "DNA", "RC", "Unknown"]
TARGET_COLORS = {"UHRF1": "#b2182b", "DNMT1": "#ef8a62", "DNMT3B": "#fddbc7", "EZH2": "#2166ac",
                 "KDM1A": "#67a9cf", "METTL3": "#1b7837", "ATRX": "#762a83"}


def bh(p):
    p = np.asarray(p, float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if ok.any():
        out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out


def te_info(feature: pd.Series) -> pd.DataFrame:
    parts = feature.astype(str).str.split(":")
    df = pd.DataFrame({"feature": feature, "is_te": parts.str.len() >= 3})
    df["family"] = parts.str[0]
    df["te_class"] = parts.str[2].where(df["is_te"]).str.rstrip("?")
    df["strict"] = df["is_te"] & ~df["te_class"].isin(NON_TE_CLASSES | {"Satellite", "Low_complexity"})
    return df


def read_counts(paths: list[Path]) -> pd.DataFrame:
    cols = {}
    for f in paths:
        t = pd.read_csv(f, sep="\t", index_col=0)
        cols[f.name.removesuffix(".cntTable")] = t.iloc[:, 0]
    m = pd.DataFrame(cols).fillna(0)
    m.index = m.index.astype(str).str.strip('"')
    return m


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--te-de", required=True)
    ap.add_argument("--counts", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--gene-names", default=None)
    ap.add_argument("--fish-05c", default=None)
    ap.add_argument("--out", required=True)
    ap.add_argument("--padj", type=float, default=0.05)
    ap.add_argument("--min-lfc", type=float, default=0.5)
    ap.add_argument("--min-base-mean", type=float, default=10)
    ap.add_argument("--min-share-change", type=float, default=5.0,
                    help="minimum |%% change| in TE share for a direction call (with FDR < 0.05)")
    args = ap.parse_args()
    out = Path(args.out)
    (out / "tables").mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(20261001)

    man = pd.read_csv(args.manifest, sep="\t")
    studies = man.drop_duplicates("study_id").set_index("study_id")
    samples = pd.read_csv(args.samples, sep="\t").set_index("sample_id")
    counts = read_counts(sorted(Path(args.counts).glob("*.cntTable")))
    info_all = te_info(pd.Series(counts.index, index=counts.index))
    strict_rows = info_all.index[info_all["strict"]]
    gene_rows = info_all.index[~info_all["is_te"]]
    names = None
    if args.gene_names:
        names = pd.read_csv(args.gene_names, sep="\t", header=None, names=["gene_id", "gene_name"]).set_index(
            "gene_id")["gene_name"]

    srows, crows, frows, grows = [], [], [], []
    for f in sorted(Path(args.te_de).glob("*.tsv")):
        sid = f.stem
        if sid not in studies.index:
            continue
        meta = studies.loc[sid]
        de = pd.read_csv(f, sep="\t")
        de["feature"] = de["feature"].astype(str).str.strip('"')
        de = de.merge(te_info(de["feature"]), on="feature")
        te = de[de["strict"] & (de["baseMean"] >= args.min_base_mean)].dropna(subset=["log2FoldChange"])
        sig = te[(te["padj"] < args.padj) & (te["log2FoldChange"].abs() >= args.min_lfc)]
        n_up, n_down = int((sig["log2FoldChange"] > 0).sum()), int((sig["log2FoldChange"] < 0).sum())
        bt = stats.binomtest(n_up, n_up + n_down) if n_up + n_down else None
        lfc = te["log2FoldChange"].to_numpy()
        boots = [np.median(rng.choice(lfc, len(lfc))) for _ in range(2000)] if len(lfc) else [np.nan]
        # global TE share of counted reads
        runs = [r for r in samples.index[samples["parental_phenotype"] == sid] if r in counts.columns]
        role = samples.loc[runs, "treatment"]
        tot = counts[runs].sum()
        share = counts.loc[strict_rows, runs].sum() / tot
        lg = np.log(share / (1 - share))
        a, b = lg[role == "perturbed"], lg[role == "control"]
        tt = stats.ttest_ind(a, b, equal_var=False)
        diff = a.mean() - b.mean()
        se = np.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b))
        dfw = se ** 4 / ((a.var(ddof=1) / len(a)) ** 2 / (len(a) - 1) + (b.var(ddof=1) / len(b)) ** 2 / (len(b) - 1))
        q = stats.t.ppf(0.975, dfw)
        inv = lambda x: 1 / (1 + np.exp(-x))  # noqa: E731
        base = b.mean()
        srows.append({
            "study_id": sid, "gse": meta["gse"], "target": meta["target_symbol"],
            "perturbation": meta["perturbation"], "sample_type": meta["sample_type"],
            "n_control": int((role == "control").sum()), "n_perturbed": int((role == "perturbed").sum()),
            "n_te_tested": len(te), "n_sig_up": n_up, "n_sig_down": n_down,
            "frac_up_sig": n_up / (n_up + n_down) if n_up + n_down else np.nan,
            "frac_up_ci95_low": bt.proportion_ci().low if bt else np.nan,
            "frac_up_ci95_high": bt.proportion_ci().high if bt else np.nan,
            "p_binom_up_vs_down": bt.pvalue if bt else np.nan,
            "median_lfc": float(np.median(lfc)) if len(lfc) else np.nan,
            "median_lfc_ci95_low": float(np.quantile(boots, 0.025)),
            "median_lfc_ci95_high": float(np.quantile(boots, 0.975)),
            "p_wilcoxon_median": float(stats.wilcoxon(lfc).pvalue) if len(lfc) > 10 else np.nan,
            "te_share_control_pct": 100 * inv(base), "te_share_perturbed_pct": 100 * inv(a.mean()),
            "te_share_rel_change_pct": 100 * (inv(base + diff) / inv(base) - 1),
            "te_share_rel_ci95_low": 100 * (inv(base + diff - q * se) / inv(base) - 1),
            "te_share_rel_ci95_high": 100 * (inv(base + diff + q * se) / inv(base) - 1),
            "te_share_logit_diff": diff, "p_te_share": float(tt.pvalue)})
        for cls, d in te.groupby("te_class"):
            if len(d) >= 5:
                crows.append({"study_id": sid, "target": meta["target_symbol"], "te_class": cls, "n": len(d),
                              "median_lfc": float(d["log2FoldChange"].median()),
                              "frac_up_sig": float((d.loc[d.padj < args.padj, "log2FoldChange"] > 0).mean())
                              if (d.padj < args.padj).any() else np.nan,
                              "p_wilcoxon": float(stats.wilcoxon(d["log2FoldChange"]).pvalue) if len(d) > 10 else np.nan})
        frows.append(te.assign(study_id=sid, target=meta["target_symbol"])[
            ["study_id", "target", "feature", "family", "te_class", "baseMean", "log2FoldChange", "padj"]])
        if names is not None:
            gde = de[~de["is_te"]].copy()
            gde["gene_name"] = gde["feature"].map(names)
            for g in ZF_GENE.get(meta["target_symbol"], []):
                r = gde[gde["gene_name"] == g]
                if len(r):
                    grows.append({"study_id": sid, "target": meta["target_symbol"], "gene": g,
                                  "perturbation": meta["perturbation"], "baseMean": float(r.baseMean.iat[0]),
                                  "log2FoldChange": float(r.log2FoldChange.iat[0]), "padj": float(r.padj.iat[0])})

    S = pd.DataFrame(srows)
    S["fdr_te_share"] = bh(S["p_te_share"])
    S["fdr_binom"] = bh(S["p_binom_up_vs_down"])
    C = pd.DataFrame(crows)
    C["fdr"] = bh(C["p_wilcoxon"])
    F = pd.concat(frows, ignore_index=True) if frows else pd.DataFrame()
    G = pd.DataFrame(grows)
    S["observed_direction"] = np.where((S["fdr_te_share"] < 0.05)
                                       & (S["te_share_rel_change_pct"].abs() >= args.min_share_change),
                                       np.sign(S["te_share_logit_diff"]), 0)

    # Expectation (silencer loss derepresses TEs) and fish ISP predictions
    S["expected_direction"] = np.where(S["perturbation"].isin(["KO", "KD"]), 1, 0)
    if args.fish_05c:
        for mode in ("plain", "rmTEpc1"):
            pf = Path(args.fish_05c) / f"05c_{mode}_per_regulator.tsv"
            if not pf.exists():
                continue
            fr = pd.read_csv(pf, sep="\t")
            preds, pps = [], []
            for r in S.itertuples():
                arm = "overexpress" if r.perturbation == "OE" else "delete"
                m = fr[(fr["perturb_type"] == arm) & fr["target"].isin(FISH_REG.get(r.target, []))]
                preds.append(float(m["frac_up_sig"].mean()) if len(m) else np.nan)
                pps.append(float(m["p_emp_hits"].min()) if len(m) else np.nan)
            S[f"fish_{mode}_pred_frac_up"] = preds
            S[f"fish_{mode}_p_emp_vs_random"] = pps
            S[f"fish_{mode}_pred_direction"] = np.sign(np.array(preds) - 0.5)
    conc = []
    obs = S[S["observed_direction"] != 0]
    for col, label in [("expected_direction", "silencer loss -> TE up (KO/KD only)")] + \
            [(f"fish_{m}_pred_direction", f"fish bulk ISP ({m})") for m in ("plain", "rmTEpc1")
             if f"fish_{m}_pred_direction" in S]:
        d = obs[obs[col] != 0].dropna(subset=[col])
        k = int((d[col] == d["observed_direction"]).sum())
        bt = stats.binomtest(k, len(d), 0.5) if len(d) else None
        ci = bt.proportion_ci(confidence_level=0.95, method="exact") if bt else None  # Clopper-Pearson
        conc.append({"prediction": label, "n_studies_with_significant_shift": len(d), "n_agree": k,
                     "frac_agree": k / len(d) if len(d) else np.nan,
                     "frac_agree_ci95_lo": ci.low if ci else np.nan,
                     "frac_agree_ci95_hi": ci.high if ci else np.nan,
                     "p_binom_vs_0.5": bt.pvalue if bt else np.nan})
    K = pd.DataFrame(conc)
    if len(K):
        pv = K["p_binom_vs_0.5"]
        K["fdr_bh"] = np.nan
        K.loc[pv.notna(), "fdr_bh"] = multipletests(pv.dropna(), method="fdr_bh")[1]

    S.to_csv(out / "tables/08_study_te_response.tsv", sep="\t", index=False)
    C.to_csv(out / "tables/08_class_te_response.tsv", sep="\t", index=False)
    F.to_csv(out / "tables/08_family_lfc.tsv.gz", sep="\t", index=False)
    G.to_csv(out / "tables/08_perturbed_gene_check.tsv", sep="\t", index=False)
    K.to_csv(out / "tables/08_prediction_concordance.tsv", sep="\t", index=False)
    figures(S, C, K, out / "figures")
    print(S[["study_id", "n_te_tested", "n_sig_up", "n_sig_down", "te_share_rel_change_pct", "fdr_te_share"]]
          .round(3).to_string(index=False))
    print(K.to_string(index=False))


# Studies whose TE result is flagged in the report: unique-mapping % differs by
# genotype (|gap| >= 5 points and |r| >= 0.7; 08c_qc_by_arm.tsv, DC-18).
FLAGGED = {"ZF_uhrf1_liver_2023": "†", "ZF_uhrf1_larva_2020": "†", "ZF_mettl3_MO_EC_2017": "†"}
# Figure 2 stars only cells whose median is large enough to read (review item):
# Wilcoxon over hundreds of shrunken log2FC values flags trivial shifts.
STAR_MIN_ABS_MEDIAN = 0.1


def _ratio_axis(ax, axis="x", lim=(0.1, 20)):
    """Log-scaled ratio axis (perturbed/control TE share) labelled as % change,
    so one very wide CI does not flatten every other study."""
    from matplotlib.ticker import FixedLocator, FuncFormatter
    ticks = [0.125, 0.25, 0.5, 1, 2, 4, 16]
    fmt = FuncFormatter(lambda v, _: f"{(v - 1) * 100:+.0f}%" if v != 1 else "0")
    if axis == "x":
        ax.set_xscale("log"); ax.set_xlim(*lim)
        ax.xaxis.set_major_locator(FixedLocator(ticks)); ax.xaxis.set_major_formatter(fmt)
        ax.xaxis.set_minor_locator(FixedLocator([]))
        ax.tick_params(axis="x", labelsize=6)
    else:
        ax.set_yscale("log"); ax.set_ylim(*lim)
        ax.yaxis.set_major_locator(FixedLocator(ticks)); ax.yaxis.set_major_formatter(fmt)
        ax.yaxis.set_minor_locator(FixedLocator([]))
        ax.tick_params(axis="y", labelsize=6)


def _ratio(pct):
    return np.maximum(1 + np.asarray(pct, dtype=float) / 100, 1e-3)


def figures(S, C, K, figs: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})
    S = S.sort_values(["target", "study_id"]).reset_index(drop=True)
    yy = np.arange(len(S))

    # Figure 1: TE share change per study (forest) + significant families up/down
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(8.6, 0.32 * len(S) + 1.6), sharey=True,
                                 gridspec_kw={"width_ratios": [1.3, 1]})
    for k, r in S.iterrows():
        col = TARGET_COLORS.get(r.target, "0.4")
        mk = {"KO": "o", "KD": "s", "OE": "^"}.get(r.perturbation, "o")
        x, lo, hi = _ratio([r.te_share_rel_change_pct, r.te_share_rel_ci95_low, r.te_share_rel_ci95_high])
        a1.errorbar(x, k, xerr=[[x - lo], [hi - x]], fmt=mk, color=col, ms=5, capsize=2, lw=1)
        a1.text(1.01, k, f"p={r.p_te_share:.2g}" + ("*" if r.fdr_te_share < 0.05 else ""),
                transform=a1.get_yaxis_transform(), fontsize=6, va="center")
        a2.barh(k, r.n_sig_up, color="#b2182b", height=0.6)
        a2.barh(k, -r.n_sig_down, color="#2166ac", height=0.6)
        if np.isfinite(r.p_binom_up_vs_down):
            a2.text(1.01, k, f"p={r.p_binom_up_vs_down:.2g}" + ("*" if r.fdr_binom < 0.05 else ""),
                    transform=a2.get_yaxis_transform(), fontsize=6, va="center")
    a1.axvline(1, color="k", lw=0.6)
    _ratio_axis(a1, "x")
    a2.axvline(0, color="k", lw=0.6)
    a1.set_yticks(yy, [f"{r.study_id.replace('ZF_', '')} ({r.perturbation}){FLAGGED.get(r.study_id, '')}"
                       for r in S.itertuples()], fontsize=6.5)
    a1.set_ylim(len(S) - 0.5, -0.5)
    a1.set_xlabel("TE share, % change vs control\n(95% CI; log ratio scale, CIs clipped)")
    a2.set_xlabel("TE families significant\ndown (blue) | up (red)")
    a1.set_title("Global TE load (Welch t on logit share)", fontsize=8)
    a2.set_title("Family-level DE (padj < 0.05, |log2FC| ≥ 0.5); binomial p", fontsize=8)
    fig.suptitle("Zebrafish TE response to loss or gain of TE regulators (* BH-FDR < 0.05 across studies)",
                 fontsize=8.5)
    fig.text(0.01, 0.005, "† Unique-mapping % differs by genotype (|gap| ≥ 5 points, |r| ≥ 0.7); TE direction "
             "not separable from library quality. See report.",
             fontsize=6)
    fig.tight_layout(rect=(0, 0.02, 0.96, 1))
    fig.subplots_adjust(wspace=0.35)
    fig.savefig(figs / "08_fig1_study_te_response.png", dpi=300)
    fig.savefig(figs / "08_fig1_study_te_response.pdf")
    plt.close(fig)

    # Figure 2: class-level median log2FC heatmap
    if len(C):
        classes = [c for c in CLASS_ORDER if c in set(C.te_class)]
        M = np.full((len(S), len(classes)), np.nan)
        star = np.full(M.shape, "", dtype=object)
        for i, r in S.iterrows():
            for j, c in enumerate(classes):
                d = C[(C.study_id == r.study_id) & (C.te_class == c)]
                if len(d):
                    M[i, j] = d.median_lfc.iat[0]
                    star[i, j] = "*" if d.fdr.iat[0] < 0.05 and abs(M[i, j]) >= STAR_MIN_ABS_MEDIAN else ""
        lim = max(0.05, np.nanquantile(np.abs(M), 0.95))
        fig, ax = plt.subplots(figsize=(4.2, 0.3 * len(S) + 1.4))
        im = ax.imshow(M, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                if star[i, j]:
                    ax.text(j, i, "*", ha="center", va="center", fontsize=8,
                            color="w" if abs(M[i, j]) > 0.6 * lim else "k")
        ax.set_xticks(range(len(classes)), classes)
        ax.set_yticks(yy, [r.study_id.replace("ZF_", "") + FLAGGED.get(r.study_id, "") for r in S.itertuples()],
                      fontsize=6.5)
        fig.colorbar(im, ax=ax, fraction=0.05, extend="both").set_label(
            f"median log2FC of families (colour clipped at ±{lim:.2f})")
        ax.set_title("TE class response (* Wilcoxon BH-FDR < 0.05 and |median| ≥ 0.1)", fontsize=8)
        fig.tight_layout()
        fig.savefig(figs / "08_fig2_class_response.png", dpi=300)
        fig.savefig(figs / "08_fig2_class_response.pdf")
        plt.close(fig)

    # Figure 3: fish ISP prediction vs zebrafish observation
    cols = [c for c in ("fish_plain_pred_frac_up", "fish_rmTEpc1_pred_frac_up") if c in S]
    if cols:
        fig, axes = plt.subplots(1, len(cols), figsize=(3.6 * len(cols) + 0.4, 3.4), squeeze=False, sharey=True)
        for ax, c in zip(axes[0], cols):
            d = S.dropna(subset=[c])
            for r in d.itertuples():
                y, lo, hi = _ratio([r.te_share_rel_change_pct, r.te_share_rel_ci95_low, r.te_share_rel_ci95_high])
                ax.errorbar(getattr(r, c), y, yerr=[[y - lo], [hi - y]],
                            fmt="o", color=TARGET_COLORS.get(r.target, "0.4"), ms=4, capsize=1.5, lw=0.8)
            ax.axhline(1, color="k", lw=0.6)
            _ratio_axis(ax, "y")
            ax.axvline(0.5, color="k", lw=0.6, ls=":")
            if len(d) > 2:
                rho, p = stats.spearmanr(d[c], d["te_share_rel_change_pct"])
                ax.set_title(f"{c.split('_')[1]}: Spearman ρ={rho:.2f}, p={p:.2g} (n={len(d)})", fontsize=7.5)
            ax.set_xlabel("fish ISP: fraction of significant\nTE families predicted up")
        axes[0][0].set_ylabel("zebrafish: TE share change\n(95% CI, log ratio scale)")
        handles = [plt.Line2D([], [], marker="o", ls="", color=v, label=k) for k, v in TARGET_COLORS.items()]
        axes[0][-1].legend(handles=handles, fontsize=6, frameon=False, loc="center left", bbox_to_anchor=(1.01, 0.5))
        fig.suptitle("Do fish bulk-ISP predictions match zebrafish knockouts?", fontsize=8.5)
        fig.tight_layout()
        fig.savefig(figs / "08_fig3_prediction_vs_observation.png", dpi=300)
        fig.savefig(figs / "08_fig3_prediction_vs_observation.pdf")
        plt.close(fig)


if __name__ == "__main__":
    main()
