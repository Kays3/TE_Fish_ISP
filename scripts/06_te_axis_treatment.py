#!/usr/bin/env python3
"""06 — Does a treatment shift the global TE-load axis of a bulk cohort?

Built for the A. polyacanthus ocean-acidification cohort (162 libraries,
brain/gill/liver x control/developmental/intergenerational x sensitive/
tolerant parental phenotype), where 05_bulk_linear_isp.py found one axis
dominating within-stratum TE variation (TE PC1 ~29%, tracking TE read
fraction). Question: is that axis also moved by CO2 exposure?

Global TE measures, per library:
  axis_score       projection onto TE PC1. Loadings are estimated from the
                   WITHIN tissue x treatment PCA (as in the ISP run), so
                   treatment differences cannot shape the axis being
                   tested; the projection uses tissue-centred data, so
                   between-treatment differences remain visible. Oriented
                   so that positive = more TE.
  te_frac_strict   logit(strict-TE counts / all counted reads)
  te_frac_all      logit(all repeat-class counts / all counted reads)

Model (frozen OA protocol, section 5A): within each tissue x parental
phenotype stratum, measure ~ C(family) + C(treatment) by OLS, with
contrasts developmental-vs-control and intergenerational-vs-control.
Family is the parental pair and is partly confounded with treatment, so
the intergenerational contrast rests mainly on family 3 (effective n ~3.6
per group). A pooled per-tissue model, ~ C(phenotype:family) +
C(treatment), is reported as well. BH-FDR within each model set.
Sensitivity: the same models + library QC covariates; the pooled model
without source-study outliers (if a `source_study_outlier` column exists).

Cross-check: is the OA paper's significant TE differential expression
aligned with this axis? Per tissue x phenotype x contrast, the fraction
of significant strict-TE features moving in the direction of their PC1
loading.

Outputs (tables_dir): 06_te_axis_scores.tsv, 06_te_axis_treatment_tests.tsv,
06_te_axis_qc_by_treatment.tsv, 06_de_alignment.tsv, 06_de_alignment_with_shift.tsv,
06_te_read_pct.tsv (TE share of reads with 95% CI), 06_phenotype_tests.tsv
(tolerant vs sensitive within tissue: mixed model with family random intercept,
and OLS on family means); figures_dir: 06_phenotype_contrasts.{png,pdf}, figures_dir:
06_te_axis_by_treatment.{png,pdf}

Usage:
  python scripts/06_te_axis_treatment.py --config runs/<run>.yaml \
      [--qc runs/qc.tsv] [--de significant-te-results.tsv]
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from utils.common import NON_TE_CLASSES, load_config, provenance, setup_logging, write_json  # noqa: E402

_spec = importlib.util.spec_from_file_location("bulk05", HERE / "05_bulk_linear_isp.py")
bulk05 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bulk05)

CONTRASTS = {"developmental_vs_control": "developmental",
             "intergenerational_vs_control": "intergenerational"}


def logit(p: pd.Series) -> pd.Series:
    p = p.clip(1e-6, 1 - 1e-6)
    return np.log(p / (1 - p))


def axis_scores(counts, samples, genes, tes, bc, log) -> tuple[pd.DataFrame, pd.Series, float]:
    # (a) loadings from the within tissue x treatment standardized TE block
    within = bulk05.normalize(counts, samples, genes, tes, dict(bc, stratify_by=["tissue", "treatment"]), log)
    te_cols = [c for c in within.columns if c in tes]
    T = within[te_cols].to_numpy()
    _, sv, vt = np.linalg.svd(T - T.mean(0), full_matrices=False)
    load = vt[0]
    if np.median(load) < 0:
        load = -load
    var1 = float(sv[0] ** 2 / (sv ** 2).sum())
    # (b) project tissue-centred/-scaled data (treatment differences retained)
    tissue = bulk05.normalize(counts, samples, genes, tes, dict(bc, stratify_by=["tissue"]), log)
    score = tissue[te_cols].to_numpy() @ load
    score = (score - score.mean()) / score.std()
    out = pd.DataFrame({"axis_score": score}, index=tissue.index)
    return out, pd.Series(load, index=te_cols), var1


def fit_contrasts(df: pd.DataFrame, measure: str, formula_rhs: str, label: dict) -> list[dict]:
    d = df.copy()
    d["treatment"] = pd.Categorical(d["treatment"], ["control", "developmental", "intergenerational"])
    fit = smf.ols(f"{measure} ~ {formula_rhs}", data=d).fit()
    rows = []
    for name, level in CONTRASTS.items():
        term = f"C(treatment)[T.{level}]"
        if term not in fit.params:
            continue
        ci = fit.conf_int().loc[term]
        resid_sd = float(np.sqrt(fit.scale))
        rows.append({**label, "measure": measure, "contrast": name,
                     "estimate": float(fit.params[term]), "se": float(fit.bse[term]),
                     "ci95_low": float(ci[0]), "ci95_high": float(ci[1]),
                     "p": float(fit.pvalues[term]), "df_resid": int(fit.df_resid),
                     "std_effect": float(fit.params[term]) / resid_sd,
                     "n": int(fit.nobs)})
    return rows


TREATS = ["control", "developmental", "intergenerational"]


def phenotype_contrasts(d: pd.DataFrame, measure: str, method: str) -> list[dict]:
    """Tolerant - sensitive within each treatment, averaged over treatments,
    and phenotype x treatment interactions (difference in CO2 effect,
    tolerant - sensitive). method: 'mixed' = random intercept per family line
    on fish-level data; 'family_means' = OLS on family x treatment means."""
    import warnings
    d = d.copy()
    d["phenotype"] = pd.Categorical(d["parental_phenotype"], ["sensitive", "tolerant"])
    d["treatment"] = pd.Categorical(d["treatment"], TREATS)
    formula = f"{measure} ~ C(phenotype) * C(treatment)"
    info = {}
    if method == "mixed":
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fit = smf.mixedlm(formula, d, groups=d["fam_id"]).fit(reml=True)
        fam_var = float(fit.cov_re.iloc[0, 0])
        info = {"converged": bool(fit.converged), "family_var": fam_var,
                "boundary_fit": bool(fam_var < 1e-6 or not fit.converged),
                "n_units": int(len(d)), "n_families": int(d["fam_id"].nunique())}
        if info["boundary_fit"]:
            # Family variance estimated at zero: the mixed model reduces to
            # fish-level OLS, which is fitted explicitly (stable SEs, t-based).
            fit = smf.ols(formula, d).fit()
            beta, cov = fit.params, fit.cov_params()
            dist, df = stats.t, fit.df_resid
            info["method_used"] = "ols_boundary"
        else:
            beta, cov = fit.fe_params, fit.cov_params().loc[fit.fe_params.index, fit.fe_params.index]
            dist, df = stats.norm, None
            info["method_used"] = "mixedlm_reml"
    else:
        fm = d.groupby(["fam_id", "phenotype", "treatment"], observed=True)[measure].mean().reset_index()
        fit = smf.ols(formula, fm).fit()
        beta, cov = fit.params, fit.cov_params()
        dist, df = stats.t, fit.df_resid
        info = {"converged": True, "family_var": np.nan, "boundary_fit": False,
                "method_used": "ols_family_means", "n_units": int(len(fm)),
                "n_families": int(fm["fam_id"].nunique())}
    names = list(beta.index)
    P = "C(phenotype)[T.tolerant]"

    def L(terms: dict) -> np.ndarray:
        v = np.zeros(len(names))
        for t, w in terms.items():
            v[names.index(t)] = w
        return v

    def inter(level):
        return f"{P}:C(treatment)[T.{level}]"
    contrasts = {
        "tolerant_minus_sensitive@control": L({P: 1}),
        "tolerant_minus_sensitive@developmental": L({P: 1, inter("developmental"): 1}),
        "tolerant_minus_sensitive@intergenerational": L({P: 1, inter("intergenerational"): 1}),
        "tolerant_minus_sensitive@mean_of_treatments": L({P: 1, inter("developmental"): 1 / 3,
                                                            inter("intergenerational"): 1 / 3}),
        "interaction_developmental_effect_tol_minus_sens": L({inter("developmental"): 1}),
        "interaction_intergenerational_effect_tol_minus_sens": L({inter("intergenerational"): 1}),
    }
    rows = []
    for name, v in contrasts.items():
        est = float(v @ beta.to_numpy())
        se = float(np.sqrt(v @ cov.to_numpy() @ v))
        stat = est / se
        if dist is stats.norm:
            p, q = float(2 * stats.norm.sf(abs(stat))), 1.959964
        else:
            p, q = float(2 * stats.t.sf(abs(stat), df)), float(stats.t.ppf(0.975, df))
        rows.append({"measure": measure, "method": method, "contrast": name, "estimate": est, "se": se,
                     "ci95_low": est - q * se, "ci95_high": est + q * se, "p": p, **info})
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--config", default=None)
    ap.add_argument("--qc", default=None, help="TSV: sample_id + numeric QC columns")
    ap.add_argument("--de", default=None, help="OA significant-te-results.tsv for the alignment check")
    args = ap.parse_args()
    cfg = load_config(args.config)
    log = setup_logging(cfg, "06_te_axis_treatment")
    bc = dict(cfg["bulk_isp"])
    bc["counts"] = cfg.path("bulk_isp", "counts")
    tables, figs = cfg.path("paths", "tables_dir"), cfg.path("paths", "figures_dir")
    tables.mkdir(parents=True, exist_ok=True)
    figs.mkdir(parents=True, exist_ok=True)

    counts = bulk05.load_counts(bc, log)
    samples = pd.read_csv(cfg.path("bulk_isp", "samples"), sep="\t", dtype=str).set_index(bc["sample_id_col"])
    counts = counts[samples.index.intersection(counts.columns)]
    samples = samples.loc[counts.columns]
    genes, tes = bulk05.split_features(counts, NON_TE_CLASSES | set(cfg.te.get("exclude_classes", [])))
    repeats = counts.index[counts.index.str.count(":") >= 2]

    sc, loadings, var1 = axis_scores(counts, samples, genes, tes, bc, log)
    total = counts.sum()
    sc["te_frac_strict"] = logit(counts.loc[tes].sum() / total)
    sc["te_frac_all"] = logit(counts.loc[repeats].sum() / total)
    sc = sc.join(samples)
    sc["fam_id"] = sc["parental_phenotype"] + "_" + sc["family"].astype(str)
    qc_cols = []
    if args.qc:
        qc = pd.read_csv(args.qc, sep="\t").set_index("sample_id")
        qc_cols = [c for c in qc.columns if pd.api.types.is_numeric_dtype(qc[c])]
        sc = sc.join(qc[qc_cols])
        for c in qc_cols:                     # standardize within tissue for the covariate models
            sc[c] = sc.groupby("tissue")[c].transform(lambda v: (v - v.mean()) / (v.std() or 1))
    sc.to_csv(tables / "06_te_axis_scores.tsv", sep="\t")
    log.info("TE PC1 (within tissue x treatment) explains %.1f%%; corr(axis_score, te_frac_strict) = %.2f",
             100 * var1, sc["axis_score"].corr(sc["te_frac_strict"]))

    rows = []
    for measure in ("axis_score", "te_frac_strict", "te_frac_all"):
        for (tissue, pheno), d in sc.groupby(["tissue", "parental_phenotype"]):
            rows += fit_contrasts(d, measure, "C(family) + C(treatment)",
                                  {"model": "stratum_family", "tissue": tissue, "phenotype": pheno})
            if qc_cols:
                rows += fit_contrasts(d, measure, "C(family) + C(treatment) + " + " + ".join(qc_cols),
                                      {"model": "stratum_family_qc", "tissue": tissue, "phenotype": pheno})
        for tissue, d in sc.groupby("tissue"):
            rows += fit_contrasts(d, measure, "C(fam_id) + C(treatment)",
                                  {"model": "tissue_pooled", "tissue": tissue, "phenotype": "both"})
            if qc_cols:
                rows += fit_contrasts(d, measure, "C(fam_id) + C(treatment) + " + " + ".join(qc_cols),
                                      {"model": "tissue_pooled_qc", "tissue": tissue, "phenotype": "both"})
    if "source_study_outlier" in sc.columns:
        for measure in ("axis_score", "te_frac_strict", "te_frac_all"):
            for tissue, d in sc[sc["source_study_outlier"].astype(str) != "yes"].groupby("tissue"):
                rows += fit_contrasts(d, measure, "C(fam_id) + C(treatment)",
                                      {"model": "tissue_pooled_no_outliers", "tissue": tissue, "phenotype": "both"})
    res = pd.DataFrame(rows)
    res["fdr"] = np.nan
    for (_, _), idx in res.groupby(["model", "measure"]).groups.items():
        res.loc[idx, "fdr"] = multipletests(res.loc[idx, "p"], method="fdr_bh")[1]
    res.to_csv(tables / "06_te_axis_treatment_tests.tsv", sep="\t", index=False)

    # Tolerant vs sensitive within tissue (family nested in phenotype)
    prow = []
    for measure in ("axis_score", "te_frac_strict", "te_frac_all"):
        for tissue, d in sc.groupby("tissue"):
            for method in ("mixed", "family_means"):
                prow += [{"tissue": tissue, **r} for r in phenotype_contrasts(d, measure, method)]
    pres = pd.DataFrame(prow)
    pres["fdr"] = np.nan
    for _, idx in pres.groupby(["method", "measure"]).groups.items():
        pres.loc[idx, "fdr"] = multipletests(pres.loc[idx, "p"], method="fdr_bh")[1]
    pres.to_csv(tables / "06_phenotype_tests.tsv", sep="\t", index=False)
    log.info("Phenotype tests: %d contrasts; mixed-model fits converged: %s",
             len(pres), bool(pres.loc[pres["method"] == "mixed", "converged"].all()))

    # Do QC metrics themselves differ by treatment? (technical explanation check)
    qrows = []
    for c in qc_cols:
        for tissue, d in sc.groupby("tissue"):
            qrows += fit_contrasts(d.rename(columns={c: "qcv"}), "qcv", "C(fam_id) + C(treatment)",
                                   {"model": "qc_by_treatment", "tissue": tissue, "phenotype": "both", "qc": c})
    if qrows:
        qres = pd.DataFrame(qrows)
        qres["fdr"] = multipletests(qres["p"], method="fdr_bh")[1]
        qres.to_csv(tables / "06_te_axis_qc_by_treatment.tsv", sep="\t", index=False)

    # Alignment of the OA paper's significant TE DE with the axis loadings
    if args.de:
        de = pd.read_csv(args.de, sep="\t")
        de = de[(de["strict_te"].astype(str).str.upper() == "TRUE") & (de["analysis_set"] == "full")]
        de = de[de["feature_id"].isin(loadings.index)]
        de["loading"] = de["feature_id"].map(loadings)
        de["aligned"] = np.sign(de["log2FoldChange"]) == np.sign(de["loading"])
        al = de.groupby(["tissue", "parental_phenotype", "contrast"]).agg(
            n_sig=("feature_id", "size"), frac_up=("log2FoldChange", lambda x: float((x > 0).mean())),
            frac_aligned_with_axis=("aligned", "mean")).reset_index()
        from scipy.stats import binomtest
        al["p_binom_aligned"] = [binomtest(int(round(r.frac_aligned_with_axis * r.n_sig)), r.n_sig, 0.5).pvalue
                                 for r in al.itertuples()]
        al["fdr"] = multipletests(al["p_binom_aligned"], method="fdr_bh")[1]
        al.to_csv(tables / "06_de_alignment.tsv", sep="\t", index=False)
        # Same, but against the direction of that stratum's estimated axis shift
        st = res[(res["model"] == "stratum_family") & (res["measure"] == "axis_score")].set_index(
            ["tissue", "phenotype", "contrast"])["estimate"]
        rows_s = []
        for (t, ph, c), d in de.groupby(["tissue", "parental_phenotype", "contrast"]):
            if c == "intergenerational_vs_developmental":
                shift = st[(t, ph, "intergenerational_vs_control")] - st[(t, ph, "developmental_vs_control")]
            elif (t, ph, c) in st.index:
                shift = st[(t, ph, c)]
            else:
                continue
            same = int((np.sign(d["log2FoldChange"]) == np.sign(shift)).sum())
            rows_s.append({"tissue": t, "phenotype": ph, "contrast": c, "axis_shift_sd": shift,
                           "n_sig": len(d), "n_same_dir_as_shift": same, "frac": same / len(d),
                           "p_binom": binomtest(same, len(d), 0.5).pvalue})
        als = pd.DataFrame(rows_s)
        als["fdr"] = multipletests(als["p_binom"], method="fdr_bh")[1]
        als.to_csv(tables / "06_de_alignment_with_shift.tsv", sep="\t", index=False)
        tot, n = int(als["n_same_dir_as_shift"].sum()), int(als["n_sig"].sum())
        log.info("DE moving with the stratum's axis shift: %d/%d (%.0f%%), binomial p=%.3g",
                 tot, n, 100 * tot / n, binomtest(tot, n, 0.5).pvalue)
        frac_pos = float((loadings > 0).mean())
        log.info("DE alignment: %d significant strict-TE rows; %.1f%% of all loadings are positive",
                 len(de), 100 * frac_pos)

    # Absolute scale: TE share of counted reads with 95% CI. Baseline = observed
    # control mean (logit); treated = baseline + family-adjusted contrast and its
    # CI from the primary pooled model (baseline uncertainty not propagated).
    def inv(x):
        return 100 / (1 + np.exp(-x))
    pct_rows = []
    for measure in ("te_frac_strict", "te_frac_all"):
        for tissue, d in sc.groupby("tissue"):
            ctl = d.loc[d["treatment"] == "control", measure]
            m, se = ctl.mean(), ctl.std(ddof=1) / np.sqrt(len(ctl))
            tcrit = float(stats.t.ppf(0.975, len(ctl) - 1))
            pct_rows.append({"measure": measure, "tissue": tissue, "group": "control",
                             "pct": inv(m), "ci95_low": inv(m - tcrit * se), "ci95_high": inv(m + tcrit * se),
                             "relative_change_pct": 0.0, "rel_ci95_low": np.nan, "rel_ci95_high": np.nan,
                             "fdr": np.nan})
            for c in res[(res["model"] == "tissue_pooled") & (res["measure"] == measure)
                         & (res["tissue"] == tissue)].itertuples():
                pct_rows.append({"measure": measure, "tissue": tissue,
                                 "group": c.contrast.replace("_vs_control", ""),
                                 "pct": inv(m + c.estimate), "ci95_low": inv(m + c.ci95_low),
                                 "ci95_high": inv(m + c.ci95_high),
                                 "relative_change_pct": 100 * (inv(m + c.estimate) / inv(m) - 1),
                                 "rel_ci95_low": 100 * (inv(m + c.ci95_low) / inv(m) - 1),
                                 "rel_ci95_high": 100 * (inv(m + c.ci95_high) / inv(m) - 1),
                                 "fdr": c.fdr})
    pd.DataFrame(pct_rows).to_csv(tables / "06_te_read_pct.tsv", sep="\t", index=False)

    # Figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    order = ["control", "developmental", "intergenerational"]
    fig, axes = plt.subplots(2, 3, figsize=(11, 6.5), sharey="row")
    rng = np.random.default_rng(0)
    for r, measure in enumerate(["axis_score", "te_frac_strict"]):
        for c, tissue in enumerate(["brain", "gill", "liver"]):
            ax = axes[r, c]
            d = sc[sc["tissue"] == tissue]
            for k, (pheno, col, off) in enumerate([("sensitive", "#c0504d", -0.15), ("tolerant", "#4f81bd", 0.15)]):
                for i, t in enumerate(order):
                    v = d[(d["treatment"] == t) & (d["parental_phenotype"] == pheno)][measure]
                    ax.scatter(i + off + rng.uniform(-0.05, 0.05, len(v)), v, s=12, color=col, alpha=0.7,
                               label=pheno if (i == 0 and c == 0 and r == 0) else None)
                    ci = float(stats.t.ppf(0.975, len(v) - 1)) * v.std(ddof=1) / np.sqrt(len(v))
                    ax.errorbar(i + off, v.mean(), yerr=ci, color=col, lw=1.5, capsize=3, zorder=3)
                    ax.hlines(v.mean(), i + off - 0.1, i + off + 0.1, color=col, lw=2)
            # model-based contrasts (primary pooled model, family-adjusted) with BH-FDR
            tm = "axis_score" if measure == "axis_score" else "te_frac_strict"
            lines = []
            for cc in res[(res["model"] == "tissue_pooled") & (res["measure"] == tm)
                          & (res["tissue"] == tissue)].itertuples():
                lab = "dev" if cc.contrast.startswith("developmental") else "intergen"
                star = "*" if cc.fdr < 0.05 else ""
                lines.append(f"{lab}: {cc.estimate:+.2f} [{cc.ci95_low:+.2f}, {cc.ci95_high:+.2f}] "
                             f"FDR {cc.fdr:.2g}{star}")
            ax.text(0.02, 0.02, "\n".join(lines), transform=ax.transAxes, fontsize=6.5,
                    va="bottom", ha="left", family="monospace",
                    bbox=dict(facecolor="white", alpha=0.8, edgecolor="none", pad=1.5))
            ax.set_xticks(range(3), ["ctrl", "dev", "intergen"])
            if r == 0:
                ax.set_title(tissue)
            if c == 0:
                ax.set_ylabel({"axis_score": "Global TE axis (SD)",
                               "te_frac_strict": "logit TE read fraction"}[measure])
    for row in axes:                          # room below the data for the statistics text
        lo, hi = row[0].get_ylim()
        row[0].set_ylim(lo - 0.32 * (hi - lo), hi)
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.suptitle("Global TE-load axis by CO2 treatment (A. polyacanthus, 162 libraries)\n"
                 "points = libraries; bars = group mean ± 95% CI; text = family-adjusted contrast vs "
                 "control [95% CI], BH-FDR (* < 0.05)", fontsize=9)
    fig.tight_layout()
    fig.savefig(figs / "06_te_axis_by_treatment.png", dpi=200)
    fig.savefig(figs / "06_te_axis_by_treatment.pdf")
    # Forest plot: tolerant - sensitive (mixed model) per tissue x treatment
    fig2, axes2 = plt.subplots(1, 3, figsize=(14, 5.2), sharey=True)
    labels = {"tolerant_minus_sensitive@control": "ctrl", "tolerant_minus_sensitive@developmental": "dev",
              "tolerant_minus_sensitive@intergenerational": "intergen",
              "tolerant_minus_sensitive@mean_of_treatments": "mean",
              "interaction_developmental_effect_tol_minus_sens": "Δdev effect",
              "interaction_intergenerational_effect_tol_minus_sens": "Δintergen effect"}
    ylabels = []
    for panel, (ax, measure) in enumerate(zip(axes2, ("axis_score", "te_frac_strict", "te_frac_all"))):
        sub = pres[(pres["method"] == "mixed") & (pres["measure"] == measure)]
        for k, (tissue, col) in enumerate([("brain", "#7f6aa8"), ("gill", "#2a9d8f"), ("liver", "#c97b2a")]):
            for j, (cname, lab) in enumerate(labels.items()):
                r = sub[(sub["tissue"] == tissue) & (sub["contrast"] == cname)].iloc[0]
                y = k * (len(labels) + 1) + j
                ax.errorbar(r.estimate, y, xerr=[[r.estimate - r.ci95_low], [r.ci95_high - r.estimate]],
                            fmt="o", color=col, ms=4, capsize=2)
                ax.text(1.02, y, f"p={r.p:.2g}" + ("*" if r.fdr < 0.05 else ""),
                        transform=ax.get_yaxis_transform(), fontsize=6.5, va="center", ha="left",
                        clip_on=False)
                if panel == 0:
                    ylabels.append((y, f"{tissue} {lab}"))
        ax.axvline(0, color="k", lw=0.7)
        ax.set_title({"axis_score": "Global TE axis (SD)", "te_frac_strict": "logit strict-TE fraction",
                      "te_frac_all": "logit all-repeat fraction"}[measure], fontsize=9)
    axes2[0].set_yticks([y for y, _ in ylabels])
    axes2[0].set_yticklabels([t for _, t in ylabels], fontsize=7)
    axes2[0].tick_params(axis="y", labelleft=True)
    axes2[0].set_ylim(max(y for y, _ in ylabels) + 0.8, -0.8)   # top-to-bottom order
    fig2.suptitle("Tolerant − sensitive parental phenotype within tissue: estimate ± 95% CI from a mixed model "
                  "with family random intercept\n(fish-level OLS where the family variance is at its boundary); "
                  "raw p shown, * BH-FDR < 0.05. Δ = phenotype × treatment interaction", fontsize=9)
    fig2.tight_layout(rect=(0, 0, 0.97, 1))
    fig2.subplots_adjust(wspace=0.32)
    fig2.savefig(figs / "06_phenotype_contrasts.png", dpi=200)
    fig2.savefig(figs / "06_phenotype_contrasts.pdf")

    write_json(tables / "06_provenance.json", provenance(
        cfg, pc1_within_var=var1, n_samples=int(len(sc)), qc_covariates=qc_cols,
        de_table=args.de))
    log.info("Wrote 06 tables and figure")


if __name__ == "__main__":
    main()
