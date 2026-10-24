#!/usr/bin/env python3
"""04 — Validate ISP predictions against public bulk perturbation RNA-seq.

Geneformer ISP predicts what deleting / overexpressing a regulator does to
TE families. Public knockout / knockdown / overexpression RNA-seq of the
same regulators is the closest thing to ground truth. PerturbAtlas
(122,801 libraries, 7,778 studies, 13 species; doi:10.1093/nar/gkae851)
and KnockTF 2.0 (doi:10.1093/nar/gkad1016) index those studies, but
neither quantifies TEs, so the raw reads are re-counted with TEcount via
TE_Ocean_Acidification/pipelines/rnaseq/align-and-count.sh.

Stages:
  select    filter the curated study manifest to configured targets,
            species and replicate floor; write one TE-pipeline sample
            sheet per study and a target-coverage table
  baseline  observational baseline: sign of single-cell co-expression
            between each regulator and each TE family (no model at all);
            a knockout is predicted to move a family opposite to rho
  compare   score ISP (03 output) and the baseline against bulk TE
            log2 fold changes: sign concordance on bulk-significant
            pairs (binomial test vs 0.5), Spearman on all pairs, and a
            regulator-level check (does the bulk TE response go the way
            the goal-state shift says?)

Between `select` and `compare`, per study:
  bash ../TE_Ocean_Acidification/pipelines/rnaseq/align-and-count.sh \
       data/bulk/samplesheets/<study>.tsv <STAR index> <gene GTF> <TE GTF> data/bulk/tecount/<study>
  Rscript scripts/04_bulk_te_deseq.R data/bulk/tecount/<study> \
       data/bulk/samplesheets/<study>.tsv data/bulk/te_de/<study>.tsv

Interpretation rule: ISP "works" for TE regulation only where it beats the
observational baseline on the same pairs. Matching bulk data is necessary,
not sufficient (cell lines differ from the single-cell context).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils.common import load_config, provenance, setup_logging, write_json  # noqa: E402


def target_symbols(cfg) -> set[str]:
    return {g for spec in cfg["targets"].values() for g in spec["genes"]}


# ---------------------------------------------------------------------------
# select
# ---------------------------------------------------------------------------
def stage_select(cfg, log) -> None:
    bv = cfg.bulk_validation
    col = bv["columns"]
    man = pd.read_csv(cfg.path("bulk_validation", "manifest"), sep="\t", dtype=str)
    missing = [v for v in col.values() if v not in man.columns]
    if missing:
        raise KeyError(f"Manifest lacks columns {missing}; map them in bulk_validation.columns")
    man = man.rename(columns={v: k for k, v in col.items()})
    man["arm"] = man["perturbation"].map(bv["arm_of"])
    wanted = target_symbols(cfg)

    def eligible(study: pd.DataFrame) -> str:
        if study["target"].iat[0] not in wanted:
            return "target_not_configured"
        if study["arm"].isna().any():
            return "perturbation_type_unmapped"
        reps = study["role"].value_counts()
        if min(reps.get("control", 0), reps.get("perturbed", 0)) < bv["min_replicates"]:
            return "too_few_replicates"
        return "eligible"

    status = man.groupby("study_id").apply(eligible).rename("status")
    man = man.join(status, on="study_id")
    man["family_level"] = man["species"].isin(bv["species"])
    sel = man[man["status"] == "eligible"]

    tables = cfg.path("paths", "tables_dir")
    tables.mkdir(parents=True, exist_ok=True)
    man.drop_duplicates("study_id")[["study_id", "target", "perturbation", "arm", "species",
                                     "sample_type", "status", "family_level"]].to_csv(
        tables / "04_study_selection.tsv", sep="\t", index=False)

    # Which configured targets have bulk ground truth, per arm?
    cov = (sel.drop_duplicates("study_id").groupby(["target", "arm"]).agg(
        n_studies=("study_id", "nunique"),
        n_human=("family_level", "sum")).reset_index())
    cov = (pd.DataFrame({"target": sorted(wanted)}).merge(cov, how="left", on="target")
           .fillna({"n_studies": 0, "n_human": 0}))
    cov.to_csv(tables / "04_target_coverage.tsv", sep="\t", index=False)
    log.info("[select] %d eligible studies for %d / %d configured targets",
             sel["study_id"].nunique(), sel["target"].nunique(), len(wanted))

    # One sample sheet per study in the TE_Ocean_Acidification manifest format:
    # sample_id read1 read2 treatment parental_phenotype strandedness
    fq = cfg.path("bulk_validation", "fastq_dir")
    out = cfg.path("bulk_validation", "tecount_dir").parent / "samplesheets"
    out.mkdir(parents=True, exist_ok=True)
    for sid, study in sel.groupby("study_id"):
        rows = []
        for r in study.itertuples():
            paired = str(r.layout).upper() == "PAIRED"
            rows.append({
                "sample_id": r.run,
                "read1": fq / (f"{r.run}_1.fastq.gz" if paired else f"{r.run}.fastq.gz"),
                "read2": fq / f"{r.run}_2.fastq.gz" if paired else "-",
                "treatment": r.role,
                "parental_phenotype": sid,          # column reused as study label
                "strandedness": r.strandedness if r.strandedness in {"no", "forward", "reverse"} else "no",
            })
        pd.DataFrame(rows).to_csv(out / f"{sid}.tsv", sep="\t", index=False)
    write_json(tables / "04_select_provenance.json", provenance(
        cfg, n_eligible=int(sel["study_id"].nunique()), samplesheets=str(out)))


# ---------------------------------------------------------------------------
# baseline
# ---------------------------------------------------------------------------
def stage_baseline(cfg, log) -> None:
    import anndata as ad
    import scanpy as sc

    proc = cfg.path("paths", "processed_dir")
    genes = ad.read_h5ad(proc / "geneformer_input.h5ad")
    te = ad.read_h5ad(proc / "te_families.h5ad")
    te = te[:, te.var["readout"].astype(bool).to_numpy()]
    train = (genes.obs[cfg.split["key"]] == "train").to_numpy()   # never look at test donors
    genes, te = genes[train].copy(), te[train].copy()
    sym = genes.var["source_id"].astype(str)
    if cfg.vocab.mode == "extended":
        is_gene = ~genes.var["ensembl_id"].astype(str).str.startswith(cfg.vocab.te_token_prefix)
        genes, sym = genes[:, is_gene.to_numpy()].copy(), sym[is_gene]
    genes.obs["n_counts_genes"] = np.asarray(genes.X.sum(1)).ravel()
    sc.pp.normalize_total(genes, target_sum=1e4)
    sc.pp.log1p(genes)

    ct_key = cfg.input.obs_keys["cell_type"]
    rows = []
    for symbol in sorted(target_symbols(cfg)):
        hits = np.flatnonzero(sym.to_numpy() == symbol)
        if not len(hits):
            continue
        x_all = np.asarray(genes.X[:, hits[0]].todense()).ravel()
        for fam_i, fam in enumerate(te.var_names):
            y_all = np.asarray(te.X[:, fam_i]).ravel()
            rhos, weights = [], []
            # within cell type, so the baseline is not just cell-type composition
            for _, idx in genes.obs.groupby(ct_key, observed=True).indices.items():
                x, y = x_all[idx], y_all[idx]
                if (x > 0).sum() < cfg.bulk_validation["baseline_min_cells_expressing"]:
                    continue
                r = spearmanr(x, y).statistic
                if np.isfinite(r):
                    rhos.append(r)
                    weights.append(len(idx))
            if rhos:
                rho = float(np.average(rhos, weights=weights))
                rows.append({"target": symbol, "family": fam, "rho": rho, "n_cell_types": len(rhos),
                             "pred_delete": -rho, "pred_overexpress": rho})
    out = pd.DataFrame(rows)
    out.to_csv(cfg.path("paths", "tables_dir") / "04_observational_baseline.tsv", sep="\t", index=False)
    log.info("[baseline] %d target x family pairs from %d train cells", len(out), train.sum())


# ---------------------------------------------------------------------------
# compare
# ---------------------------------------------------------------------------
def te_family(feature: str) -> str:
    """TEcount names TE rows 'family:superfamily:class'."""
    parts = str(feature).split(":")
    return parts[0] if len(parts) >= 3 else feature


def load_bulk_de(cfg, log) -> pd.DataFrame:
    bv = cfg.bulk_validation
    sel = pd.read_csv(cfg.path("paths", "tables_dir") / "04_study_selection.tsv", sep="\t")
    sel = sel[sel["status"] == "eligible"].set_index("study_id")
    frames = []
    for f in sorted(cfg.path("bulk_validation", "te_de_dir").glob("*.tsv")):
        sid = f.stem
        if sid not in sel.index:
            log.warning("TE DE table %s has no eligible study row — ignored", f.name)
            continue
        de = pd.read_csv(f, sep="\t")
        de = de[de["feature"].astype(str).str.count(":") >= 2]          # TE rows only
        de = de[de["baseMean"] >= bv["min_base_mean"]]
        de["family"] = de["feature"].map(te_family)
        de["study_id"] = sid
        for k in ("target", "arm", "species", "family_level"):
            de[k] = sel.loc[sid, k]
        frames.append(de)
    if not frames:
        raise FileNotFoundError("No per-study TE DE tables — run TEcount + 04_bulk_te_deseq.R")
    de = pd.concat(frames, ignore_index=True)
    de["bulk_sig"] = (de["padj"] < bv["de_padj"]) & (de["log2FoldChange"].abs() >= bv["de_min_abs_lfc"])
    return de


def concordance(df: pd.DataFrame, pred: str, truth: str = "log2FoldChange") -> dict:
    d = df.dropna(subset=[pred, truth])
    sig = d[d["bulk_sig"]]
    agree = int((np.sign(sig[pred]) == np.sign(sig[truth])).sum())
    res = {"n_pairs": len(d), "n_bulk_sig": len(sig), "sign_agree": agree,
           "sign_concordance": agree / len(sig) if len(sig) else np.nan,
           "p_binom_vs_0.5": binomtest(agree, len(sig), 0.5, alternative="greater").pvalue if len(sig) else np.nan,
           "spearman_all": spearmanr(d[pred], d[truth]).statistic if len(d) > 2 else np.nan}
    return res


def stage_compare(cfg, log, isp_table: str | None = None) -> None:
    tables = cfg.path("paths", "tables_dir")
    de = load_bulk_de(cfg, log)

    # Pool replicate studies per (target, arm, family): mean LFC, any-significant.
    human = de[de["family_level"].astype(bool)]
    pooled = (human.groupby(["target", "arm", "family"]).agg(
        log2FoldChange=("log2FoldChange", "mean"), bulk_sig=("bulk_sig", "any"),
        n_studies=("study_id", "nunique")).reset_index())

    isp_table = Path(isp_table) if isp_table else tables / "03_regulator_te_family_effects.tsv"
    log.info("[compare] ISP predictions from %s", isp_table)
    isp = pd.read_csv(isp_table, sep="\t")
    isp = isp[isp["status"] == "tested"].rename(columns={"perturb_type": "arm", "effect_vs_null": "isp"})
    base = pd.read_csv(tables / "04_observational_baseline.tsv", sep="\t").melt(
        id_vars=["target", "family"], value_vars=["pred_delete", "pred_overexpress"],
        var_name="arm", value_name="baseline")
    base["arm"] = base["arm"].str.replace("pred_", "", regex=False)

    m = (pooled.merge(isp[["target", "arm", "family", "isp", "fdr"]], how="left",
                      on=["target", "arm", "family"])
         .merge(base, how="left", on=["target", "arm", "family"]))
    m.to_csv(tables / "04_bulk_vs_predictions.tsv", sep="\t", index=False)

    both = m.dropna(subset=["isp", "baseline"])            # same pairs for a fair contest
    rows = []
    for arm, sub in [("all", both)] + list(both.groupby("arm")):
        for name in ("isp", "baseline"):
            rows.append({"arm": arm, "predictor": name, **concordance(sub, name)})
    scores = pd.DataFrame(rows)
    scores.to_csv(tables / "04_family_level_scores.tsv", sep="\t", index=False)

    # Regulator level (all species): does the study's TE response go where
    # the goal-state shift says? Bulk side = fraction of significant TE
    # families that went up; ISP side = sign of shift toward TE_high.
    reg = (de[de["bulk_sig"]].groupby(["study_id", "target", "arm", "species"])
           .agg(frac_up=("log2FoldChange", lambda x: float((x > 0).mean())),
                n_sig_te=("family", "nunique")).reset_index())
    shift_f = tables / "03_goal_state_shift_vs_null.tsv"
    if shift_f.exists() and isp_table is None:
        sh = pd.read_csv(shift_f, sep="\t").rename(columns={"perturb_type": "arm"})
        reg = reg.merge(sh[["target", "arm", "effect_vs_null", "fdr"]], how="left", on=["target", "arm"])
        reg["agree"] = np.sign(reg["frac_up"] - 0.5) == np.sign(reg["effect_vs_null"])
    reg.to_csv(tables / "04_regulator_level.tsv", sep="\t", index=False)

    isp_all = scores[(scores["arm"] == "all") & (scores["predictor"] == "isp")].iloc[0] if len(scores) else None
    base_all = scores[(scores["arm"] == "all") & (scores["predictor"] == "baseline")].iloc[0] if len(scores) else None
    verdict = "insufficient pairs"
    if isp_all is not None and isp_all["n_bulk_sig"] >= 10:
        verdict = ("ISP beats observational baseline"
                   if isp_all["sign_concordance"] > base_all["sign_concordance"]
                   and isp_all["p_binom_vs_0.5"] < 0.05
                   else "ISP does NOT beat observational baseline")
    lines = [
        "# Bulk ground-truth validation", "",
        f"**Verdict: {verdict}**", "",
        f"Human studies with TE DE tables: {human['study_id'].nunique()}; "
        f"target x arm x family pairs scored by both predictors: {len(both)}.", "",
        scores.to_markdown(index=False, floatfmt=".3g") if len(scores) else "_no overlapping pairs_", "",
        "## Reading this", "",
        "- Sign concordance is computed only on pairs the bulk experiment calls significant; "
        "0.5 is chance.",
        "- The baseline uses no model and no perturbation: single-cell co-expression in "
        "TRAIN donors. If it matches bulk as well as ISP does, ISP adds nothing for this question.",
        "- Bulk experiments are mostly cell lines with near-total knockouts or acute "
        "knockdowns; disagreement can be context, not model failure.",
    ]
    (tables / "04_bulk_validation.md").write_text("\n".join(lines) + "\n")
    write_json(tables / "04_compare_provenance.json", provenance(cfg, verdict=verdict))
    log.info("[compare] %s", verdict)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("stage", choices=("select", "baseline", "compare"))
    ap.add_argument("--config", default=None)
    ap.add_argument("--isp-table", default=None,
                    help="predictions to score (default: 03 Geneformer output; "
                         "e.g. results/tables/05_bulk_isp_effects.tsv)")
    args = ap.parse_args()
    cfg = load_config(args.config)
    log = setup_logging(cfg, "04_bulk_validation")
    if args.stage == "compare":
        stage_compare(cfg, log, args.isp_table)
    else:
        {"select": stage_select, "baseline": stage_baseline}[args.stage](cfg, log)


if __name__ == "__main__":
    main()
