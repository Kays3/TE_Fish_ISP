# Simple in silico perturbation for bulk RNA-seq with genes and TEs

`scripts/05_bulk_linear_isp.py` is a CPU-only, CellOracle-style linear
network ISP for bulk libraries in which genes and TE families are counted
together (TEcount `.cntTable`, as produced by the TE_Ocean_Acidification
pipeline).

## Why this method

| Option | Problem for bulk + TEs |
|---|---|
| Geneformer ISP on bulk profiles | Pretrained on single cells. A bulk mixture is out of distribution. Human genes only and no TE tokens. |
| BulkFormer | No in silico perturbation support, and no TE features |
| Deep perturbation predictors (GEARS, scGPT, …) | Need Perturb-seq training data; did not beat simple baselines in benchmarks ([Ahlmann-Eltze et al. 2025](https://doi.org/10.1038/s41592-025-02772-6)) |
| Co-expression alone (`04 baseline`) | Marginal correlation. It cannot separate a regulator from its correlated neighbours, and it has no indirect effects. |
| **Linear network + propagation (this script)** | **Simple and interpretable. TEs are just more target nodes, and the method works for any species (no tokenizer, no orthology), so it applies directly to the 162-library *A. polyacanthus* cohort.** |

The recipe is CellOracle's (Kamimoto *et al.* 2023, *Nature* 614:742,
[doi:10.1038/s41586-022-05688-9](https://doi.org/10.1038/s41586-022-05688-9)):
a ridge-regression network fit on unperturbed data, a perturbation clamped
on one regulator, and the shift propagated for 3 steps. CellOracle restricts
edges with a chromatin-accessibility base GRN, which bulk RNA-seq lacks.
Here every configured regulator may connect to every TE family, and
shrinkage plus a permutation null do the regularising.

## Steps

1. **Normalize.** DESeq2 median-of-ratios size factors are computed on
   *gene* rows only, so TE multimappers do not set library size. Then
   log2(count/sf + 1), with a low-count filter.
2. **Stratify.** Every feature is centred and scaled within each stratum
   (`stratify_by`, e.g. tissue × treatment), so the network learns
   co-variation within a tissue rather than differences between tissues.
3. **Network.** A closed-form multi-output ridge maps regulators to
   {regulators, strict-TE families}, with no self-edges and one penalty
   chosen by CV. Held-out R² is recorded per TE family; families below
   `min_cv_r2` are not reported.
4. **Perturb.** The regulator is clamped at its 1st percentile (delete) or
   99th percentile (overexpress). The shift propagates through the
   regulator→regulator block for `n_steps` and then onto TEs. Step 1 (direct
   edges only) is always reported alongside.
5. **Uncertainty.** A stratified sample bootstrap gives a CI and sign
   stability. A permutation null shuffles the regulator within stratum and
   refits the network, giving an empirical p and BH-FDR.
6. **Validation.** The output has the same columns as the Geneformer
   results, so it can be scored against PerturbAtlas ground truth with:

```bash
python scripts/04_bulk_validation.py compare --isp-table results/tables/05_bulk_isp_effects.tsv
```

## Synthetic check (planted answers; 162 libraries, 3 tissues)

The synthetic data has the following planted effects:

- trim28 represses L1-1;
- setdb1a acts on L1-1 only *through* trim28;
- dnmt1 represses ERV1-2;
- yy1 activates L2-3.

L1-1 also carries a tissue effect that runs *opposite* to its within-tissue
repression, and a simple repeat is linked to trim28.

| Check | Result |
|---|---|
| All planted edges, both arms, right sign | Yes. 8/8 at FDR 0.0037 (1000 permutations); bootstrap sign stability 1.0 |
| Indirect regulator | setdb1a → L1-1 shift is 0.38 with direct edges only and 1.38 after 3 steps |
| False positives | 1 of 10 hits (trim28 → L2-3, FDR 0.042), about what a 5% FDR allows |
| Simple repeat | Dropped by the strict-TE filter |
| **Without tissue stratification** | **trim28 → L1-1 sign flips (−1.9 instead of +1.8), and every random TE family gets a large spurious effect. CV R² still looks fine (0.39), so R² does not catch this.** |
| Runtime | About 5 s on a laptop (5 regulators, 1000 permutations, 200 bootstraps) |

## Applying it to the fish cohort

- **Counts.** Point `bulk_isp.counts` at the TEcount `.cntTable` directory
  for the 162-library run on the NAS, and `samples` at a sheet with
  `sample_id`, `tissue` and `treatment`.
- **Stratification.** Use `stratify_by: [tissue, treatment]` to keep CO2
  exposure from posing as regulation. That leaves about 18 libraries per
  stratum (9 per phenotype within each). The script warns below 10.
- **Regulators.** The configured symbols match teleost paralogs
  automatically (`setdb1a`/`setdb1b`). Add
  `extra_regulator_patterns: ["^znf\\d+$"]` for the KRAB-like ZNFs that
  stood out in the TE↔gene bridge, then check `05_bulk_target_cv_r2.tsv`
  before reading any effect.
- **Scale.** n ≈ 160 supports tens of regulators, not hundreds. With p ≳ n
  the ridge still fits, but effects shrink toward zero and the permutation
  null decides.

## Applied to the fish cohort (2026-10-01): negative

On the 162-library *A. polyacanthus* cohort, the 30 regulators present in the
annotation were indistinguishable from expression-matched random genes,
both before and after removing the dominant global TE-load axis. Details,
configs and tables are in
[`runs/2026-10-01_fish_bulk_isp/RUN_RECORD.md`](runs/2026-10-01_fish_bulk_isp/RUN_RECORD.md).
**Always run `--random-regulators` controls**: without them the plain run
reported 32,142 "significant" regulator → TE effects.

## Limits

- This is a linear model of observational co-variation. A predicted effect
  is a hypothesis; the step-1 effect is a partial correlation, not a
  causal estimate.
- A regulator that barely varies across samples cannot show an effect,
  however important it is.
- Unmeasured common drivers (cell composition, stress) produce edges that
  stratification cannot remove. Check top hits against the
  PerturbAtlas/KnockTF ground truth (`04`) before treating them as
  regulation.
