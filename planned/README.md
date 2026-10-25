# Planned work (not yet run)

Nothing in this folder has been run, downloaded or analysed in this repository. It reserves the layout for the next
two lines of work so that they can reuse the fish ISP code, the manifests and the figure driver.

| Folder | Question | Status |
|---|---|---|
| `acute_te_response/` | How do TE families respond within days of CO2 or temperature exposure, and can a perturbation model say which regulators account for it? | planned, not yet run |
| `warming/` | The same question for ocean warming in the clownfish *Amphiprion ocellaris*, and for warming against the zebrafish background. | planned, not yet run |

The completed TE-expression analyses both folders start from (acute CO2 in zebrafish; warming in *A. ocellaris* with its rRNA re-check), with their scripts, tables and reports, are in the private repository [Kays3/TE_Acute_Warming](https://github.com/Kays3/TE_Acute_Warming) (`main` at `7049e41`).

Rules both inherit from the fish work (see `docs/history.md`):

1. Register the primary test, sidedness, predicted direction, alpha, minimum attainable p, the multiple-testing
   family and the status vocabulary (ISP-STD-1 section A) before any output exists, and commit it first.
2. Beat expression-matched random genes (B1), never zero alone, with at least 20 matched controls per gene (the fish run had 3 per regulator, which left per-gene status `not_estimable`); score delete and overexpress on the identical sample
   list; report donors, animals or tanks as the unit, not libraries.
3. Include positive controls (C2) and state power by simulation (C3). Without a positive control that responds,
   a negative cannot be read, and the status is `not_estimable`.
4. Inputs by path and sha256 under `manifests/`; large data stays outside git; no GPU unless a plan names it.
