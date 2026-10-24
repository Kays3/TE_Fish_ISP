# Acute TE response (planned, not yet run)

**Status: planned. No analysis has been run in this repository.**

## What exists elsewhere (private `TE_Ocean_Acidification` repository, reports of October 2026)

- An acute CO2 check in adult zebrafish brain and gill (public cohort PRJNA1102381, 26 libraries) under a
  pre-registered design: report `zf-acute-20261006`.
- The design note for using zebrafish as a background for TE responses to CO2 and warming in damselfish:
  `zf-background-te-design-20261005`.

These are TE-expression analyses, not ISP. Their results are not restated here; read them in their own reports.

## What this folder is for

Whether a regulator-level perturbation model (the CellOracle-style ISP in `scripts/05_bulk_linear_isp.py` or a
successor) can say anything about an acute response, given what the fish ISP run showed: unperturbed variation in
regulator mRNA did not separate real regulators from matched random genes (`docs/history.md`, section 1). Any plan
here must say what design avoids that failure before it runs.

## Layout (empty until the plan is registered)

```
acute_te_response/
  PREREGISTRATION.md     registered before any output; append-only amendments
  config/                run configs (paths through ${DATA_ROOT})
  manifests/             sha256 of every input
  results/               tables written by the registered scripts
  figures/               rebuilt by script from results/
```

## Open before anything runs

- The question (ISP, or TE expression only) and the data set. Data not already on disk need approval before download.
- The positive control that would show the design can detect a real effect.
- Compute: CPU only unless the plan states otherwise.
