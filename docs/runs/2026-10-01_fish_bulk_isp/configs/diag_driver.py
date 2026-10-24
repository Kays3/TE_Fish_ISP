import sys, importlib.util, numpy as np, pandas as pd, yaml
sys.path.insert(0, "scripts")
spec = importlib.util.spec_from_file_location("m", "scripts/05_bulk_linear_isp.py"); m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
from utils.common import load_config, NON_TE_CLASSES
import logging; log = logging.getLogger("d"); logging.basicConfig(level=logging.WARNING)
cfg = load_config("runs/fish_bulk_isp_20261001.yaml"); bc = dict(cfg["bulk_isp"])
counts = m.load_counts(bc, log)
samples = pd.read_csv("runs/fish_samples.tsv", sep="\t", dtype=str).set_index("sample_id")
counts = counts[samples.index.intersection(counts.columns)]; samples = samples.loc[counts.columns]
genes, tes = m.split_features(counts, NON_TE_CLASSES)
x = m.normalize(counts, samples, genes, tes, bc, log)
regs = m.pick_regulators(x.columns, cfg, bc, log)
te_cols = [c for c in x.columns if c in tes]
gene_cols = [c for c in x.columns if c in genes]
T = x[te_cols].to_numpy(); u, s, vt = np.linalg.svd(T - T.mean(0), full_matrices=False)
pc1 = u[:, 0] * s[0]; var1 = s[0]**2 / (s**2).sum()
# orient PC1 so that most TE loadings are positive
if np.median(vt[0]) < 0: pc1 = -pc1; vt = -vt
print(f"TE block: {len(te_cols)} families; PC1 explains {var1:.1%} of within-stratum TE variance; "
      f"{(vt[0] > 0).mean():.1%} of families load positively")
G = x[gene_cols].to_numpy(); ug, sg, _ = np.linalg.svd(G - G.mean(0), full_matrices=False)
print(f"Gene block PC1 explains {sg[0]**2/(sg**2).sum():.1%}; |corr(TE PC1, gene PC1)| = {abs(np.corrcoef(pc1, ug[:,0])[0,1]):.2f}")
te_frac = counts.loc[tes].sum() / counts.sum()
qc = pd.read_csv("runs/qc.tsv", sep="\t").set_index("sample_id").reindex(x.index)
def wcorr(v):  # within-stratum correlation with PC1
    v = pd.Series(np.asarray(v, float), index=x.index)
    grp = samples[["tissue", "treatment"]].astype(str).agg("|".join, axis=1)
    vz = v.groupby(grp).transform(lambda a: (a - a.mean()) / (a.std() or 1))
    return float(np.corrcoef(pc1, vz)[0, 1])
rows = {"TE read fraction": te_frac.reindex(x.index), "library size (total counts)": counts.sum().reindex(x.index)}
for c in qc.columns: rows[c] = pd.to_numeric(qc[c], errors="coerce")
print("\nWithin-stratum correlation with TE PC1:")
for k, v in rows.items(): print(f"  {k:32s} {wcorr(v): .2f}")
print("\nRegulators, correlation with TE PC1 (sorted):")
cr = pd.Series({r: np.corrcoef(pc1, x[r])[0, 1] for r in regs}).sort_values()
print(cr.round(2).to_string())
rg = pd.Series({g: np.corrcoef(pc1, x[g])[0, 1] for g in np.random.default_rng(0).choice(gene_cols, 2000, replace=False)})
print(f"\n2000 random genes: |corr| median {rg.abs().median():.2f}, 90th pct {rg.abs().quantile(.9):.2f}; regulators |corr| median {cr.abs().median():.2f}")
