# Fish ISP history, 1 to 13 October 2026

What the floor ran, what it measured and what it could not show. There was one fish in silico perturbation (ISP)
run, a bulk linear ISP on the *Acanthochromis polyacanthus* ocean-acidification (OA) cohort, and one ground-truth
study that tested it against 17 zebrafish knockout, knockdown and overexpression contrasts. Every number below is
quoted from a committed report or table in `results/`; where a claim was not checked it is marked. Times are JST.

Result wording follows ISP-STD-1 v1.1 (status vocabulary: `positive`, `negative`, `opposite_direction`,
`not_estimable`, `control_draw_sensitive_open`, `no_op_failed`, `not_run_by_design`, `stopped_not_analysed`). The
fish run is described by the criterion it did not meet (B1, matched controls), not as an absence of an effect.

## At a glance

| Study | Reading | Source |
|---|---|---|
| Fish bulk ISP, 30 TE regulators | B1 not met: the regulators were not distinguishable from 90 expression-matched random genes, plain or with the TE axis removed | `results/2026-10-01_fish_bulk_isp/RUN_RECORD.md`, `tables/05c_*` |
| Zebrafish ground truth, 17 contrasts | 0 of 16 loss-of-function studies shifted global TE share at FDR < 0.05; weakly powered | `results/2026-10-01_zebrafish_validation/REPORT.md` |
| Fish prediction against zebrafish response | Not supported: Spearman rho = -0.18 (p = 0.50, n = 17); the test crosses species, tissue and perturbation type | `tables/08_prediction_concordance.tsv`, `08f_exploratory_concordance.tsv` |

The fish run used 0 GPU-hours. The standard cited is ISP-STD-1 v1.1; the fish work predates the E3 cross-evaluation
of the criteria and was read against them afterwards.

## 1. The fish bulk ISP run

On 1 October the bulk linear ISP script (`scripts/05_bulk_linear_isp.py`, a CellOracle-style ridge network with three
propagation steps) was run on the 162-library TEcount gene and TE matrix of the *A. polyacanthus* cohort, stratified
by tissue and CO2 treatment. Thirty configured TE regulators were present in the fish annotation. TRIM28/KAP1 has no
annotated orthologue, and the KRAB-ZNF set, DUX4, SUV39H1/2 and TEX19 are absent, so the human positive controls did
not transfer. The method is closed-form linear algebra on a 162 x ~2,200 matrix; no GPU was used.

The decisive step was the random-gene control: 90 expression-matched pseudo-regulators over three seeds. Without it
the plain run would have looked like 32,142 significant regulator-to-TE effects (28% of 114,720 tests). With it, the
real regulators produced no more hits than random genes (median FDR < 0.05 families 393 against 371, Mann-Whitney
p = 0.70), and one regulator reached empirical p <= 0.05 where 1.5 were expected by chance. No regulator passed BH
correction.

The diagnostics explain why. A single global TE-load axis dominates the within-stratum TE variance (PC1 29.3%,
PC2 3.5%), 96.8% of families load with the same sign, and the axis tracks TE read fraction at r = 0.87 [0.83, 0.91].
In the plain run each regulator's predicted direction follows its correlation with that axis: deleting *tasora*
(+0.64) lowers almost every family and deleting *tp53* (-0.45) raises them. Removing the axis left 3,728 of 59,520
tests at FDR < 0.05 and still no separation from random genes (p = 0.17). The *atrx* deletion predicted LINE
families up (14 of 16 significant LINE hits), which fits ATRX's known role in repeat silencing, but its empirical p
was 0.055 and BH-FDR 0.82, and the run record does not report it as a finding.

The run record's reading is that natural variation in silencer mRNA across unperturbed animals carries no detectable
information about TE regulation beyond what any gene of similar expression carries, and that bulk linear ISP cannot
nominate TE regulators in this dataset. Under ISP-STD-1 this is a failure of the matched-control rule (B1). It does
not show that these regulators fail to act on fish TEs: transcript level is a poor proxy for silencing activity, and
the design cannot separate the two.

| Quantity | Value |
|---|---|
| Libraries | 162 (154 without outliers) |
| Regulators / random pseudo-regulators | 30 / 90 (3 seeds) |
| Plain: hits real vs random | 393 vs 371, Mann-Whitney p 0.70 |
| Plain: median abs effect | Mann-Whitney p 0.74 |
| Empirical p <= 0.05 | 1 of 30 (expected 1.5); binomial excess p 0.79 |
| Axis removed: hits | 62.5 vs 68, Mann-Whitney p 0.17 |
| TE PC1 variance | 29.3% (PC2 3.5%) |
| PC1 vs TE read fraction | r 0.87 [0.83, 0.91] |
| Outlier robustness | rho 0.985; Jaccard 0.80 |
| Review | no named reviewer in the run record |

A side result from the same cohort, outside ISP: developmental CO2 exposure lowered the TE read share in gill by
11.4% [5.1, 17.2] (FDR 0.005); brain and liver were inconclusive
(`results/2026-10-01_te_axis_treatment/REPORT.md`). What the TE axis is biologically (nascent or intronic RNA
fraction, or cell composition) was left open.

Figure: `figures/2026-10-01_fish_bulk_isp/05_fig1_real_vs_random.png` (new in this repository; plots the saved 05c
values).

## 2. The zebrafish ground truth

The study asked whether knocking out, knocking down or overexpressing a TE regulator in zebrafish moves TE
expression, and whether the direction matches two predictions: that losing a silencer derepresses TEs, and the
direction the fish ISP run gave for the matching regulator. It ran no ISP itself. The panel was 17 contrasts from 14
GEO series and 112 bulk RNA-seq runs (uhrf1, dnmt1, dnmt3, kdm1a, ezh2, mettl3, atrx, and human UHRF1
overexpression), aligned with STAR 2.7.11b and counted with TEcount 2.2.4 in multi mode against 3,562,104 danRer11 TE
loci. The tests were fixed in code at `bd7086c` of the source repository (1 Oct 19:05 JST) before the first count
table existed (about 23:10 JST).

No knockout or knockdown study changed global TE share at FDR < 0.05 (0 of 16). The one study that moved was human
UHRF1 overexpression in larval liver, +29.9% [16.6, 44.5], FDR 0.011, and the report reads it as compositional:
gene-assigned reads fell from 10.27 to 7.73 million per library (-2.55 [-4.98, -0.11], p = 0.043) while
TE-assigned reads did not change (p = 0.86). With no loss-of-function study shifted, the prespecified silencer-loss
concordance test had nothing to test.

The family-level unit, defined after the prespecified tests came back near empty and so exploratory, leans one way.
Among loss-of-function studies with a significant lean, 6 of 8 leaned up (0.75 [0.35, 0.97], p = 0.29); among DNA
methylation regulators 6 of 7 (p = 0.13); and 5 of 5 once the three studies with a genotype-linked mapping gap are
dropped (p = 0.063). None passes FDR. The direction is consistent, the evidence is weak, and two of the leans
coincide with library-quality differences.

| Quantity | Value |
|---|---|
| Panel | 17 contrasts, 14 GEO series, 112 runs |
| KO/KD TE-share shifts | 0 of 16 at FDR < 0.05 |
| hUHRF1 OE TE share | +29.9% [16.6, 44.5], FDR 0.011 |
| hUHRF1 OE gene reads | -2.55 M [-4.98, -0.11], p 0.043 |
| Family lean, KO/KD | 6 of 8 up, p 0.29 |
| Family lean, DNA methylation | 6 of 7 up, p 0.13 |
| GSE234993 unique mapping | 20.6% mutant vs 26.9% control (read 1) |
| Review | science and numbers review passed at `246c3e4` and `95031c8` (source repository), final |

Figures 1 to 3 of that report are `figures/2026-10-01_zebrafish_validation/08_fig{1,2,3}_*.png`.

## 3. Prediction against ground truth

The fish ISP directions did not track the zebrafish responses. The Spearman correlation between the predicted
fraction of TE families up and the observed TE-share change was -0.18 (p = 0.50, n = 17) for the plain run and -0.11
(p = 0.67) with the TE axis removed. On family-level leans the plain predictions agreed in 4 of 9 studies (p = 1.0).
The one prespecified test with any study in it found 0 of 1 agreeing (Clopper-Pearson 0 to 0.975).

This is consistent with the run record, since the fish calls were already no better than random genes. It is also a
weak test. The step crosses species (*A. polyacanthus* to *Danio rerio*), tissue, and perturbation type
(environmental CO2 to knockout), and several zebrafish studies share one fish prediction (four uhrf1 studies, five
dnmt1 studies), so the number of independent regulator-level tests is smaller than the study counts suggest. The
comparison therefore only weakly disconfirms the predictions.

## 4. What these studies cannot show

Most zebrafish arms have 2 to 5 libraries, so a non-shift in one study is weak evidence against a response, and
global TE share is a compositional readout that moves when gene reads change. Counting is unstranded, so antisense
transcription through a TE counts as TE expression, and multi-mode counting spreads reads across near-identical
young subfamilies. Morpholinos and two alleles (uhrf1 hi272, dnmt1 s904) do not lower their own mRNA, so the
perturbed-gene check cannot confirm the perturbation there. The positive controls, dnmt1 and uhrf1 loss, did not
respond at the TE-share level, which under ISP-STD-1 C2 limits what any negative can mean. No attainability or power check was run before the study (A2, C3). Clutch blocking was
added after an interim look (DC-19) and created one of the six upward leans. In GSE234993 the unique-mapping gap
between mutants and controls persists on read 1 alone, and that study's TE result does not survive adjustment for
mapping rate (post hoc, DC-18).

On the fish side, the cohort is unperturbed animals, the regulator signal is confounded with a global TE axis whose
biology is unknown, and removing that axis is known from synthetic data to add false positives as well as remove
confounding.

## 5. Later readings, dropped items and what was not verified

- The 2 October synthesis and its 12 October refresh carried the zebrafish result into the cross-programme reading:
  measuring TE derepression after regulator perturbation is harder than assumed, and the insensitivity sits in the
  summary statistic more than in the data, since the dnmt1 knockout liver shows 112 TE families up against 26 down.
- On 5 October the raw data folder (1,992 files, 834.6 GB) moved to a GPU host and the alignment-host copy was
  deleted. Reads, alignments and reference exist in one copy; counts and QC are also in the archive that
  `manifests/inputs.sha256.tsv` references.
- The 162-library cohort is unpublished. The fish and OA analyses were dropped from the 13 October lab-talk list for
  that reason.
- The `atrx` LINE lean (14 of 16) in the fish run is recorded but was never reported as a finding; it did not beat
  the random-gene control.
- No fish ISP was rerun after 2 October. The run record's follow-ups (a regulator-specific design from public
  perturbation data or a perturbation experiment in fish cells, and reruns removing two or three TE components with
  matched random controls) were not started.
- The 13 October talk list marked the zebrafish "6 of 8 / rho -0.18" pair as not verified, because the gate had not
  opened the source. Both numbers are in `results/2026-10-01_zebrafish_validation/REPORT.md` as quoted above.
  This history was compiled from committed material only; nothing was recomputed.

## 6. Relation to the E3 criteria work

The fish work is not part of E3. E3 tested the ISP-STD-1 criteria on Geneformer and human lung data (steps 1 to 11)
and closed on 24 October 2026 with its internal final report. The fish studies predate it and were read against
ISP-STD-1 v1.1 after the fact: the zebrafish report maps itself to the standard (met: tests fixed before output,
studies as the unit, opposite directions not counted; not met: attainability and power, A2 and C3), and the fish ISP
run is the floor's clearest case of why B1 exists. Without the matched random genes, 32,142 significant effects
would have been reported. E3's finding that sensitivity was never tested for any criterion applies here too:
neither fish study had a positive control that responded, so neither can say how large a real effect it would have
detected.

## Provenance of this page

Compiled on 24 October 2026 from the fish ISP history page of the same date, which was itself drawn from committed
material of the source repository (`docs/runs/2026-10-01_fish_bulk_isp/RUN_RECORD.md`, commits `526cc05` to
`2ad71cd`; `docs/runs/2026-10-01_zebrafish_validation/REPORT.md`, final `246c3e4`, last touched `74cc729`).
See `docs/source_map.md` for where each file of this repository came from.
