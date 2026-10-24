#!/usr/bin/env python3
"""05 — Simple in silico perturbation for BULK RNA-seq with genes AND TEs.

A CellOracle-style linear network ISP (Kamimoto et al. 2023, Nature,
doi:10.1038/s41586-022-05688-9), adapted to bulk samples with TE families as
additional target nodes. No GPU, no foundation model: numpy + scikit-learn,
minutes on a few hundred libraries.

  1. Counts: TEcount .cntTable files (genes and TEs in one table, as produced
     by TE_Ocean_Acidification/pipelines/rnaseq/align-and-count.sh) or one
     feature x sample TSV. TE rows are "family:superfamily:class"; non-TE
     repeat classes are dropped (strict_te).
  2. Normalize: median-of-ratios size factors computed on GENE rows only,
     log2(count / sf + 1), low-count filter, then each feature is centred
     (and scaled) WITHIN stratum (e.g. tissue), so the network learns
     within-tissue co-variation rather than tissue identity.
  3. Network: multi-output ridge, regulators -> {regulators, genes, TEs}, one
     shared penalty chosen by cross-validation. Self-edges are removed.
     Held-out R2 per target is recorded; poorly predicted targets are not
     reported.
  4. Perturb: a regulator is clamped at its stratum 1st percentile (delete)
     or 99th percentile (overexpress); the shift propagates through the
     regulator -> regulator block for `n_steps` (CellOracle uses 3), then
     onto genes and TEs. Step 1 = direct effects only.
  5. Uncertainty: sample bootstrap (stratified) -> CI and sign stability.
     Null: the regulator's values are permuted within stratum and the
     network refit -> empirical p per regulator x TE family; BH-FDR.
  6. Outputs use the same columns as 03 (target, perturb_type, family,
     effect_vs_null, fdr, status), so 04_bulk_validation.py compare
     --isp-table results/tables/05_bulk_isp_effects.tsv scores this method
     against PerturbAtlas ground truth and the co-expression baseline.

Interpretation: this is a linear model of observational co-variation. A
predicted effect is a hypothesis about regulation, with the same limits as
the co-expression baseline plus conditioning on the other regulators. It
cannot see a regulator that does not vary across the samples.

Usage:
  python scripts/05_bulk_linear_isp.py [--config config/config.yaml]
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import KFold
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils.common import NON_TE_CLASSES, load_config, provenance, setup_logging, write_json  # noqa: E402


# ---------------------------------------------------------------------------
# 1-2. Load and normalize
# ---------------------------------------------------------------------------
def load_counts(bc, log) -> pd.DataFrame:
    src = Path(bc["counts"])
    if src.is_dir():
        cols = {}
        for f in sorted(src.rglob("*.cntTable")):
            t = pd.read_csv(f, sep="\t", index_col=0)
            cols[f.name.removesuffix(".cntTable")] = t.iloc[:, 0]
        if not cols:
            raise FileNotFoundError(f"No .cntTable files under {src}")
        counts = pd.DataFrame(cols).fillna(0)
    else:
        counts = pd.read_csv(src, sep="\t", index_col=0)
    # TEcount quotes gene IDs and NCBI-derived GTFs prefix them ("gene-setdb1a");
    # TE rows stay "family:superfamily:class".
    prefix = bc.get("gene_id_prefix", "gene-")
    names = counts.index.astype(str).str.strip('"')
    names = names.where(names.str.count(":") >= 2,
                        names.str.replace(rf"^{re.escape(prefix)}", "", regex=True) if prefix else names)
    counts.index = names
    counts = counts.groupby(level=0).sum()
    log.info("Counts: %d features x %d samples", *counts.shape)
    return counts


def split_features(counts: pd.DataFrame, exclude_classes: set[str]) -> tuple[pd.Index, pd.Index]:
    parts = counts.index.to_series().str.split(":")
    is_te = parts.str.len() >= 3
    te_class = parts.str[2].where(is_te)
    strict = is_te & ~te_class.isin(exclude_classes)
    return counts.index[~is_te], counts.index[strict]


def strata_of(samples: pd.DataFrame, bc) -> pd.Series:
    """One stratum per combination of the stratify_by columns (e.g. tissue x treatment)."""
    cols = list(bc.get("stratify_by") or [])
    if not cols:
        return pd.Series("all", index=samples.index)
    return samples[cols].astype(str).agg("|".join, axis=1)


def normalize(counts, samples, genes, tes, bc, log) -> pd.DataFrame:
    g = counts.loc[genes]
    pos = g[(g > 0).all(axis=1)]
    if len(pos) < 100:
        raise ValueError("Too few genes expressed in every sample for size factors")
    logpos = np.log(pos)
    sf = np.exp((logpos.sub(logpos.mean(axis=1), axis=0)).median(axis=0))   # DESeq2 median-of-ratios
    keep = (counts.ge(bc["min_count"]).sum(axis=1) >= bc["min_samples"])
    feats = counts.index[keep & counts.index.isin(genes.union(tes))]
    x = np.log2(counts.loc[feats].div(sf, axis=1) + 1).T                     # samples x features
    grp = strata_of(samples.loc[x.index], bc)
    x = x.groupby(grp).transform(lambda v: v - v.mean())
    if bc.get("scale", True):
        sd = x.groupby(grp).transform("std").replace(0, np.nan)
        x = (x / sd).fillna(0)
    log.info("Normalized: %d samples x %d features (%d TE families)", x.shape[0], x.shape[1],
             x.columns.isin(tes).sum())
    return x


def pick_regulators(columns: pd.Index, cfg, bc, log) -> list[str]:
    """Configured target symbols, matched case-insensitively with an optional
    teleost paralog suffix a/b and NCBI copy suffix (setdb1a, dnmt3bb.1), plus
    extra regexes. Restricting the suffix to a/b keeps mycn/mych from posing
    as MYC."""
    symbols = {g for spec in cfg["targets"].values() for g in spec["genes"]}
    pats = [re.compile(rf"^{re.escape(s)}[ab]?(\.\d+)?$", re.I) for s in symbols]
    pats += [re.compile(p, re.I) for p in bc.get("extra_regulator_patterns", [])]
    regs = [c for c in columns if ":" not in c and any(p.match(c) for p in pats)]
    log.info("Regulators found in data: %d (%s)", len(regs), ", ".join(regs[:20]))
    if not regs:
        raise ValueError("No configured regulator present in the count matrix")
    found = {s for s in symbols if any(re.match(rf"^{re.escape(s)}[ab]?(\.\d+)?$", r, re.I) for r in regs)}
    log.info("Configured targets absent from this annotation: %s", ", ".join(sorted(symbols - found)))
    return regs


def matched_random_regulators(counts: pd.DataFrame, x: pd.DataFrame, regs: list[str],
                              genes: pd.Index, seed: int, n_bins: int = 20) -> list[str]:
    """Negative control: as many non-regulator genes as there are regulators,
    matched on mean log expression (quantile bins). If these pseudo-regulators
    'regulate' TEs as strongly as the real ones, the real effects mostly track
    a sample-level common driver rather than specific regulation."""
    pool = [c for c in x.columns if c in genes and c not in regs]
    level = np.log2(counts.loc[pool + regs].mean(axis=1) + 1)
    bins = pd.qcut(level.rank(method="first"), n_bins, labels=False)
    rng = np.random.default_rng(seed)
    chosen = []
    for r in regs:
        cand = [g for g in pool if bins[g] == bins[r] and g not in chosen]
        chosen.append(str(rng.choice(cand)))
    return chosen


def remove_te_axes(x: pd.DataFrame, te_cols: list[str], k: int, log) -> pd.DataFrame:
    """Regress the top-k principal components of the (already within-stratum
    standardized) TE block out of EVERY feature, regulators included.

    Bulk TE data are often dominated by one global TE-load axis (e.g. nascent /
    intronic RNA fraction) that most genes also track; without this, a
    regulator's predicted effect is mostly its correlation with that axis.
    After removal the model asks about family-specific effects only, so a
    regulator that truly moves global TE load is no longer detectable."""
    T = x[te_cols].to_numpy()
    u, sv, _ = np.linalg.svd(T - T.mean(0), full_matrices=False)
    F = u[:, :k]
    explained = (sv[:k] ** 2).sum() / (sv ** 2).sum()
    X = x.to_numpy()
    X = X - F @ (F.T @ X)
    log.info("Removed %d TE principal component(s) (%.1f%% of TE-block variance) from all features",
             k, 100 * explained)
    return pd.DataFrame(X, index=x.index, columns=x.columns)


# ---------------------------------------------------------------------------
# 3. Network
# ---------------------------------------------------------------------------
def fit_network(x: pd.DataFrame, regs: list[str], targets: list[str], alpha: float) -> np.ndarray:
    """W (regulators x targets) by closed-form multi-output ridge, fitting each
    regulator-as-target without itself (no self-loops)."""
    X = x[regs].to_numpy()
    Y = x[targets].to_numpy()
    xtx = X.T @ X
    W = np.linalg.solve(xtx + alpha * np.eye(len(regs)), X.T @ Y)
    for j, t in enumerate(targets):            # regulator targets: drop the self predictor
        if t in regs:
            i = regs.index(t)
            others = [k for k in range(len(regs)) if k != i]
            Xo = X[:, others]
            W[:, j] = 0
            W[others, j] = np.linalg.solve(Xo.T @ Xo + alpha * np.eye(len(others)), Xo.T @ Y[:, j])
    return W


def choose_alpha(x, regs, targets, bc, seed) -> tuple[float, pd.Series]:
    ridge = RidgeCV(alphas=bc["alphas"], cv=KFold(bc["cv_folds"], shuffle=True, random_state=seed))
    ridge.fit(x[regs].to_numpy(), x[targets].to_numpy())
    alpha = float(ridge.alpha_)
    # held-out R2 per target at that alpha
    pred = np.zeros((len(x), len(targets)))
    for tr, te in KFold(bc["cv_folds"], shuffle=True, random_state=seed).split(x):
        W = fit_network(x.iloc[tr], regs, targets, alpha)
        pred[te] = x.iloc[te][regs].to_numpy() @ W
    y = x[targets].to_numpy()
    r2 = 1 - ((y - pred) ** 2).sum(0) / ((y - y.mean(0)) ** 2).sum(0)
    return alpha, pd.Series(r2, index=targets)


# ---------------------------------------------------------------------------
# 4. Perturbation by propagation
# ---------------------------------------------------------------------------
def propagate(W: np.ndarray, regs: list[str], targets: list[str], reg: str, delta: float,
              n_steps: int) -> np.ndarray:
    """CellOracle-style: v0 = delta on the perturbed regulator; each step
    v_reg <- v_reg @ W_RR with the perturbed regulator clamped; final shift on
    every target = v_reg @ W. Returns the shift per target."""
    r_idx = [targets.index(r) for r in regs]
    W_rr = W[:, r_idx]
    i = regs.index(reg)
    v = np.zeros(len(regs))
    v[i] = delta
    for _ in range(n_steps - 1):
        v = v @ W_rr
        v[i] = delta
    shift = v @ W
    shift[targets.index(reg)] = delta
    return shift


def perturbation_delta(x: pd.DataFrame, reg: str, arm: str, q: dict) -> float:
    """Mean distance from current level to the clamp level (stratum-scaled units)."""
    v = x[reg].to_numpy()
    level = np.quantile(v, q["delete"] if arm == "delete" else q["overexpress"])
    return float(np.mean(level - v))


# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--config", default=None)
    ap.add_argument("--random-regulators", type=int, default=None, metavar="SEED",
                    help="negative control: replace regulators with expression-matched random "
                         "genes; tables are written with prefix 05_random<SEED>_")
    args = ap.parse_args()
    cfg = load_config(args.config)
    log = setup_logging(cfg, "05_bulk_linear_isp")
    bc = dict(cfg["bulk_isp"])
    bc["counts"] = cfg.path("bulk_isp", "counts")
    rng = np.random.default_rng(cfg.project["seed"])
    tables = cfg.path("paths", "tables_dir")
    tables.mkdir(parents=True, exist_ok=True)

    counts = load_counts(bc, log)
    samples = pd.read_csv(cfg.path("bulk_isp", "samples"), sep="\t", dtype=str).set_index(bc["sample_id_col"])
    if bc.get("exclude_where"):
        col, val = bc["exclude_where"]["column"], str(bc["exclude_where"]["value"])
        drop = samples.index[samples[col].astype(str) == val]
        log.info("Excluding %d samples where %s == %s", len(drop), col, val)
        samples = samples.drop(drop)
    shared = counts.columns.intersection(samples.index)
    if len(shared) < 20:
        raise ValueError(f"Only {len(shared)} samples shared between counts and sample sheet")
    counts, samples = counts[shared], samples.loc[shared]
    genes, tes = split_features(counts, NON_TE_CLASSES | set(cfg.te.get("exclude_classes", [])))
    x = normalize(counts, samples, genes, tes, bc, log)
    regs = pick_regulators(x.columns, cfg, bc, log)
    k_pcs = int(bc.get("remove_te_pcs", 0) or 0)
    if k_pcs:
        x = remove_te_axes(x, [c for c in x.columns if c in tes], k_pcs, log)
    tag = "05_" if not k_pcs else f"05_rmTEpc{k_pcs}_"
    if args.random_regulators is not None:
        real = regs
        regs = matched_random_regulators(counts, x, regs, genes, args.random_regulators)
        tag = tag.replace("05_", f"05_random{args.random_regulators}_", 1)
        log.warning("NEGATIVE CONTROL: %d expression-matched random genes replace %s",
                    len(regs), ", ".join(real))
    te_cols = [c for c in x.columns if c in tes]
    # Non-regulator genes never feed back into propagation, so only regulators
    # and TE families need to be network targets (keeps bootstrap/permutation cheap).
    targets = regs + te_cols

    alpha, r2 = choose_alpha(x, regs, targets, bc, cfg.project["seed"])
    W = fit_network(x, regs, targets, alpha)
    rho = np.max(np.abs(np.linalg.eigvals(W[:, :len(regs)])))
    log.info("Ridge alpha=%g; spectral radius of regulator block=%.3f; TE families with "
             "CV R2 >= %.2f: %d / %d", alpha, rho, bc["min_cv_r2"],
             int((r2[te_cols] >= bc["min_cv_r2"]).sum()), len(te_cols))
    if rho >= 1:
        log.warning("Regulator block spectral radius >= 1: multi-step propagation amplifies; "
                    "prefer n_steps=1 or a larger alpha")

    p_floor = 1 / (bc["n_permutations"] + 1)
    log.info("Permutation p-value floor = %.4f (n_permutations=%d); raise n_permutations if "
             "FDR cannot reach %.2f across the tests run", p_floor, bc["n_permutations"],
             cfg.analysis["fdr"])
    strata = strata_of(samples.loc[x.index], bc)
    small = strata.value_counts()
    if (small < 10).any():
        log.warning("Strata with < 10 samples: %s — centring them removes most signal",
                    small[small < 10].to_dict())
    arms = bc["perturb_types"]
    te_idx = [targets.index(t) for t in te_cols]

    def effects(xm: pd.DataFrame, Wm: np.ndarray) -> dict:
        out = {}
        for reg in regs:
            for arm in arms:
                d = perturbation_delta(xm, reg, arm, bc["clamp_quantiles"])
                for k in sorted({1, bc["n_steps"]}):
                    out[(reg, arm, k)] = propagate(Wm, regs, targets, reg, d, k)[te_idx]
        return out

    obs = effects(x, W)

    # Bootstrap (resample samples within stratum) -> CI and sign stability
    boot = {key: [] for key in obs}
    for _ in range(bc["n_bootstrap"]):
        idx = np.concatenate([rng.choice(np.flatnonzero(strata.to_numpy() == s),
                                         (strata == s).sum()) for s in strata.unique()])
        xb = x.iloc[idx]
        for key, val in effects(xb, fit_network(xb, regs, targets, alpha)).items():
            boot[key].append(val)

    # Permutation null: regulator's values shuffled within stratum, network refit
    null = {key: [] for key in obs}
    for reg in regs:
        for _ in range(bc["n_permutations"]):
            xp = x.copy()
            for s in strata.unique():
                m = (strata == s).to_numpy()
                xp.loc[m, reg] = rng.permutation(xp.loc[m, reg].to_numpy())
            Wp = fit_network(xp, regs, targets, alpha)
            for arm in arms:
                d = perturbation_delta(xp, reg, arm, bc["clamp_quantiles"])
                for k in sorted({1, bc["n_steps"]}):
                    null[(reg, arm, k)].append(propagate(Wp, regs, targets, reg, d, k)[te_idx])

    rows = []
    for (reg, arm, k), eff in obs.items():
        b = np.array(boot[(reg, arm, k)])
        n = np.array(null[(reg, arm, k)])
        for j, fam in enumerate(te_cols):
            p = (1 + np.sum(np.abs(n[:, j]) >= abs(eff[j]))) / (1 + len(n))
            rows.append({
                "target": reg, "perturb_type": arm, "n_steps": k,
                "family": fam.split(":")[0], "te_feature": fam,
                "te_class": fam.split(":")[2],
                "effect": eff[j],
                "ci95_low": np.quantile(b[:, j], 0.025), "ci95_high": np.quantile(b[:, j], 0.975),
                "sign_stability": float(np.mean(np.sign(b[:, j]) == np.sign(eff[j]))),
                "null_mean": float(n[:, j].mean()),
                "effect_vs_null": eff[j] - float(n[:, j].mean()),
                "p_perm": p, "family_cv_r2": float(r2[fam]),
            })
    res = pd.DataFrame(rows)
    res["status"] = np.where(res["family_cv_r2"] >= bc["min_cv_r2"], "tested", "poorly_predicted")
    res["fdr"] = np.nan
    for k in res["n_steps"].unique():
        m = (res["status"] == "tested") & (res["n_steps"] == k)
        if m.any():
            res.loc[m, "fdr"] = multipletests(res.loc[m, "p_perm"], method="fdr_bh")[1]
    res.to_csv(tables / f"{tag}bulk_isp_all_steps.tsv", sep="\t", index=False)
    final = res[res["n_steps"] == bc["n_steps"]]
    final.to_csv(tables / f"{tag}bulk_isp_effects.tsv", sep="\t", index=False)

    # Network edges (direct regulator -> TE weights) for inspection
    edges = pd.DataFrame(W[:, te_idx], index=regs, columns=te_cols).stack().rename("weight").reset_index()
    edges.columns = ["regulator", "te_feature", "weight"]
    edges.to_csv(tables / f"{tag}bulk_network_reg_to_te.tsv", sep="\t", index=False)
    r2.rename("cv_r2").to_csv(tables / f"{tag}bulk_target_cv_r2.tsv", sep="\t")

    sig = final[(final["status"] == "tested") & (final["fdr"] < cfg.analysis["fdr"])]
    write_json(tables / f"{tag}bulk_isp_provenance.json", provenance(
        cfg, alpha=alpha, spectral_radius=float(rho), n_samples=int(len(x)),
        n_regulators=len(regs), regulators=regs, n_te_families=len(te_cols),
        n_significant=int(len(sig)), counts=str(bc["counts"])))
    log.info("Done: %d regulator x arm x TE tests at %d steps, %d with FDR < %.2f",
             int((final["status"] == "tested").sum()), bc["n_steps"], len(sig), cfg.analysis["fdr"])


if __name__ == "__main__":
    main()
