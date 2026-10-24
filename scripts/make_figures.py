#!/usr/bin/env python3
"""Rebuild the ten committed figures of the fish ISP work from the saved result tables in docs/runs/*/tables/.

No counts, no NAS, no GPU: tables in, figures out. Takes seconds.

  python scripts/make_figures.py [--out rebuilt] [--check]

Writes <out>/<run>/<name>.png and .pdf (06: 2 figures, 07: 5, 08: 3). With --check, compares each rebuilt PNG with the
PNG committed in docs/runs/<run>/figures/ (sha256) and exits 1 on any difference. Byte identity holds only with the
pinned environment (environment/README.md): matplotlib draws a different layout in other minor versions.
Figure code is imported from the original analysis scripts where they expose it (07 make_figures, 08 figures)
and from fig06_plots.py, a verbatim copy of the inline plotting block of 06, where they did not.
The 05 run drew no figure and none is invented here.
"""
from __future__ import annotations

import argparse
import importlib.util
import os
from pathlib import Path

os.environ.setdefault("SOURCE_DATE_EPOCH", "1791000000")  # fixed PDF creation date (2026-10-03 UTC)

import pandas as pd
import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RES = ROOT / "docs" / "runs"  # committed tables live next to the run records, as in the source repository


def _load(name: str, file: str):
    spec = importlib.util.spec_from_file_location(name, HERE / file)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def fig06(out: Path) -> None:
    t = RES / "2026-10-01_te_axis_treatment" / "tables"
    sc = pd.read_csv(t / "06_te_axis_scores.tsv", sep="\t", index_col=0)
    res = pd.read_csv(t / "06_te_axis_treatment_tests.tsv", sep="\t")
    pres = pd.read_csv(t / "06_phenotype_tests.tsv", sep="\t")
    d = out / "2026-10-01_te_axis_treatment"; d.mkdir(parents=True, exist_ok=True)
    _load("fig06_plots", "fig06_plots.py").plot_06(sc, res, pres, d)


def fig07(out: Path) -> None:
    t = RES / "2026-10-01_gene_axis" / "tables"
    rd = lambda n: pd.read_csv(t / f"{n}.tsv", sep="\t")
    gsets = yaml.safe_load(open(ROOT / "config" / "gene_sets_apoly.yaml"))
    d = out / "2026-10-01_gene_axis"; d.mkdir(parents=True, exist_ok=True)
    _load("s07", "07_treatment_gene_axis.py").make_figures(
        rd("07_gene_treatment_effects"), rd("07_set_treatment_shift"), rd("07_gene_axis_association"),
        rd("07_gene_mediation"), rd("07_set_mediation"), gsets, d)


def fig08(out: Path) -> None:
    t = RES / "2026-10-01_zebrafish_validation" / "tables"
    S = pd.read_csv(t / "08_study_te_response.tsv", sep="\t")
    C = pd.read_csv(t / "08_class_te_response.tsv", sep="\t")
    K = pd.read_csv(t / "08_prediction_concordance.tsv", sep="\t")
    d = out / "2026-10-01_zebrafish_validation"; d.mkdir(parents=True, exist_ok=True)
    _load("s08", "08_zebrafish_validation.py").figures(S, C, K, d)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=str(ROOT / "rebuilt"))
    ap.add_argument("--check", action="store_true", help="compare rebuilt PNGs with the committed PNGs")
    args = ap.parse_args()
    out = Path(args.out)
    for f in (fig06, fig07, fig08):
        f(out)
        print("done", f.__name__)
    if args.check:
        import hashlib
        h = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
        bad = 0
        for p in sorted(out.glob("*/*.png")):
            ref = RES / p.parent.name / "figures" / p.name
            ok = ref.exists() and h(p) == h(ref)
            bad += not ok
            print(("identical  " if ok else "DIFFERENT  ") + f"{p.parent.name}/{p.name}")
        print(f"{bad} of {len(list(out.glob('*/*.png')))} PNGs differ from the committed figures")
        raise SystemExit(1 if bad else 0)


if __name__ == "__main__":
    main()
