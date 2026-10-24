# Environment

| File | What it pins | Used by |
|---|---|---|
| `python-bulk.lock.txt` | The Python 3.11.15 environment that produced the 2026-10-01 results and figures (`uv pip freeze`: numpy 2.4.6, pandas 3.0.6, scipy 1.17.1, scikit-learn 1.9.1, statsmodels 0.15.0, matplotlib 3.11.2, pyyaml 6.0.3, tabulate 0.10.0). | scripts 04 to 07, 08 to 08f, `make_figures.py` |
| `deseq-frozen.packages.txt` | Package list of the frozen DESeq2 environment (R 4.4.3, DESeq2 1.46.0, apeglm 1.28.0; 191 conda packages, names with version and build). | `scripts/04_bulk_te_deseq.R`, `04b_...R`, `pipelines/zebrafish_validation/zf_deseq_*.sh` |

Alignment and counting used STAR 2.7.11b, TEcount 2.2.4 and fastp from the `TE_Ocean_Acidification` pipeline
(`align-and-count.sh` is copied into `pipelines/zebrafish_validation/`); that environment is not rebuilt here.

Create the Python environment:

```bash
uv venv --python 3.11 .venv-bulk
uv pip install --python .venv-bulk/bin/python -r environment/python-bulk.lock.txt
PYTHON=.venv-bulk/bin/python scripts/run_all.sh
```

**Figures are byte-reproducible only with this lock.** Matplotlib changes tick, font and layout details between
minor versions: the same tables drawn with matplotlib 3.10.8 give a different PNG (same data, shifted layout), while
matplotlib 3.11.2 gives the committed bytes. PDFs carry a fixed creation date (`SOURCE_DATE_EPOCH`) so they repeat
too. Installing the lock is a download; ask before doing it on a machine that does not already have the packages.
