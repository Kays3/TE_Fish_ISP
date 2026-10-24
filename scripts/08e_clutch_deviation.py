#!/usr/bin/env python3
"""Show the clutch-blocking deviation (DC-19): ~ block + role (primary) vs
~ role for the three clutch-paired studies, with the 08 family rules.

Usage: 08e_clutch_deviation.py IN_DIR OUT_TSV   (IN_DIR holds te_de/ and te_de_noblock/)
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd
from statsmodels.stats.multitest import multipletests

spec = importlib.util.spec_from_file_location("s08b", Path(__file__).with_name("08b_gse234993_sensitivity.py"))
s08b = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s08b)

STUDIES = ["ZF_uhrf1_liver_2020", "ZF_uhrf1_liver_2023", "ZF_hUHRF1_OE_liver_2023"]


def main(inp: Path, out: Path) -> None:
    rows = [{"study_id": s, "design": d, **s08b.summarise(inp / sub / f"{s}.tsv")}
            for s in STUDIES for d, sub in (("~ block + role (primary)", "te_de"), ("~ role", "te_de_noblock"))]
    T = pd.DataFrame(rows)
    T["fdr_bh_binom"] = multipletests(T["p_binom_up_vs_down"].fillna(1), method="fdr_bh")[1]
    T.to_csv(out, sep="\t", index=False)
    print(T.round(4).to_string(index=False))


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
