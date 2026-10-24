# Bulk perturbation data: PerturbAtlas and related approaches

**Question.** Can public bulk RNA-seq of genetic perturbations, PerturbAtlas
first, strengthen or replace the single-cell Geneformer ISP for TE
regulation?

**Short answer.** Use it to **validate** the ISP, not as model input.
PerturbAtlas is the best index of real knockout, knockdown and overexpression
experiments for our target regulators. It does not quantify TEs, so the
useful object is the raw reads behind it, re-counted with the TEcount
pipeline we already run for TE_Ocean_Acidification. Those real experiments
then score the ISP predictions against a baseline that uses no model.
This is implemented in `scripts/04_bulk_validation.py`.

Status: design plus implemented scaffold. No study has been downloaded or
scored yet.

## 1. What PerturbAtlas is

Zhang *et al.*, *Nucleic Acids Research* 53:D1112 (2025),
[doi:10.1093/nar/gkae851](https://doi.org/10.1093/nar/gkae851).

| | |
|---|---|
| Scale | 122,801 RNA-seq libraries, 7,778 studies, 18,612 perturbation-vs-control comparisons, 13 species (human, mouse, rat, zebrafish, chicken, pig, fly, worm, yeast, *E. coli*, *Arabidopsis*, rice, tomato) |
| Sources | GEO, ArrayExpress, SRA, ENCODE, curated through February 2023. Studies need at least one biological replicate per arm. |
| Perturbations | Knockout, knockdown (si/shRNA, CRISPRi), knock-in, mutation, overexpression |
| Processing | Trimmomatic → STAR 2.7.10b → **Salmon** 1.10.1 (gene and transcript TPM/counts) → DESeq2 (DEGs), DRIMSeq (transcript usage), rMATS (splicing) |
| Access | Web browser (forward search by target/tissue, reverse search by gene), GSEA against all signatures, and bulk downloads. The site at `perturbatlas.kratoss.site` returned HTTP 503 on 2026-10-01; the page lists Google Drive, Baidu and TeraBox mirrors. Pipeline: [ygidtu/perturb_processing](https://github.com/ygidtu/perturb_processing). Harmonizome mirrors the DEG signatures as gene sets ([human](https://maayanlab.cloud/Harmonizome/dataset/PerturbAtlas+Signatures+of+Differentially+Expressed+Genes+for+Gene+Perturbations), [mouse](https://maayanlab.cloud/Harmonizome/dataset/PerturbAtlas+Signatures+of+Differentially+Expressed+Genes+for+Mouse+Gene+Perturbations)). |

**The key limitation for this project: there are no TE counts.** Salmon
quantification against a gene/transcript annotation discards or misassigns
repeat-derived reads, so PerturbAtlas tables cannot show whether a TRIM28
knockout derepressed L1HS. Its metadata (perturbed gene, perturbation type,
tissue, run accessions) is exactly what is needed to find the experiments
that can show it.

## 2. Options considered

| # | Approach | What it gives | Verdict |
|---|---|---|---|
| A | **PerturbAtlas → raw reads → TEcount → per-study TE DE**, then score ISP and a no-model baseline | Real TE responses to real perturbations of our targets, and a fair test of whether ISP adds anything | **Adopt.** Implemented as `04_bulk_validation.py` |
| B | PerturbAtlas DEG tables as gene-level ground truth | Checks whether ISP moves the right *genes* (e.g. KRAB-ZNFs, ISGs after TRIM28 loss) with no re-alignment | **Adopt as a cheap second check.** Readout genes in `config.readout_genes` can be compared directly |
| C | PerturbAtlas reverse search and GSEA: submit the TE_high signature *genes* and ask which perturbations produce it | Candidate upstream regulators we did not think of, to add to `targets` | **Use for target discovery.** Gene-level only; output is a hypothesis list |
| D | KnockTF 2.0 (Feng *et al.*, *NAR* 52:D183, 2024, [doi:10.1093/nar/gkad1016](https://doi.org/10.1093/nar/gkad1016)): 1,468 curated knockdown/knockout datasets of 612 TFs and 172 co-factors | Overlaps PerturbAtlas for TF and co-factor targets (YY1, MYC, TRIM28 as a co-factor) | **Supplementary study source** for the same manifest |
| E | Tokenize bulk samples as pseudo-cells and run Geneformer ISP on them | Model predictions in the same context as the ground truth | **Do not adopt.** Geneformer was pretrained on single cells. A bulk profile is a mixture with a different rank structure, so the shifts would be out-of-distribution and uninterpretable. |
| F | BulkFormer: a bulk-transcriptome foundation model, ~150M parameters, ~20k protein-coding genes, >500k human bulk profiles ([bioRxiv 2025.06.11.659222](https://www.biorxiv.org/content/10.1101/2025.06.11.659222v1); [code](https://github.com/KangBoming/BulkFormer), MIT) | A bulk-native embedding | **Not for ISP.** Its README states it does not support in silico gene perturbation, and the authors removed their perturbation benchmark because of a modality mismatch. No TE features. Could serve as a bulk embedding for option B in future. |
| G | Regulon-activity inference on bulk (decoupleR / VIPER-style, using TF→target priors) | Upstream TF activity behind a bulk TE response | **Optional.** A useful cross-check on TF targets; uninformative for chromatin silencers such as HUSH and DNMT1, which have no TF-style regulons |
| H | LINCS / CMap L1000 knockdown signatures | Thousands of shRNA/CRISPR knockdowns | **Not useful.** About 1,000 measured landmark genes and no TE signal |

## 3. Why validation needs a baseline

Ahlmann-Eltze, Huber & Anders, *Nature Methods* 22:1657 (2025),
[doi:10.1038/s41592-025-02772-6](https://doi.org/10.1038/s41592-025-02772-6),
compared five single-cell foundation models and two other deep models against
deliberately simple baselines for predicting perturbation effects. None beat
the baselines. Agreement with bulk data alone therefore does not show that
Geneformer ISP adds anything.

`04_bulk_validation.py baseline` builds the simplest competitor: Spearman
co-expression between each regulator and each TE family, within cell type,
in **training donors only**, with no model and no perturbation. A knockout
is predicted to move a family opposite to ρ. `compare` scores ISP and this
baseline **on the same target × arm × family pairs**, and the verdict line
says which one wins.

## 4. Pipeline

```
PerturbAtlas / KnockTF metadata ─► data/bulk/perturbation_studies.tsv  (curated by hand)
        │  04 select: configured targets, species, ≥2 replicates/arm
        ▼
data/bulk/samplesheets/<study>.tsv   (TE_Ocean_Acidification manifest format)
        │  align-and-count.sh: fastp → STAR → TEcount --mode multi  (human TE GTF)
        ▼
data/bulk/tecount/<study>/*.cntTable
        │  04_bulk_te_deseq.R: DESeq2 ~ role, apeglm shrinkage
        ▼
data/bulk/te_de/<study>.tsv
        │  04 baseline (single-cell train donors)    03 outputs (ISP)
        ▼                                            ▼
                04 compare ─► 04_bulk_validation.md, 04_family_level_scores.tsv,
                              04_regulator_level.tsv, 04_bulk_vs_predictions.tsv
```

Two levels of scoring:

- **Family level (human studies only).** For each target × arm × TE family,
  the mean bulk log2FC across studies is compared with the ISP effect (target
  minus matched null) and with the baseline. Metrics are sign concordance on
  bulk-significant pairs (binomial test against 0.5) and Spearman correlation
  across all pairs.
- **Regulator level (any species).** For each study, the fraction of
  significant TE families that went up is compared with the sign of the ISP
  shift toward TE_high. This is how mouse and zebrafish studies contribute:
  their TE families (IAP, MERVL; teleost-specific families) have no human
  counterpart, but "TRIM28 loss derepresses TEs" is a statement that can be
  checked in any species.

## 5. Starting the manifest

1. Search PerturbAtlas (forward search) or KnockTF for each symbol in
   `config.targets`, for human first and then mouse. Mouse ESC knockouts of
   *Trim28*, *Setdb1* and *Dnmt1* are the classic TE-derepression
   experiments and are likely well represented.
2. Record one row per run, with the columns named in
   `bulk_validation.columns`. Prefer studies whose control is matched
   (same cell line and batch) and with at least two replicates per arm.
3. Fetch the FASTQ files by run accession (e.g. `fasterq-dump`) into
   `data/bulk/fastq/`. This is the expensive step: budget for the STAR index
   and tens of GB per study.
4. Run `04 select`, then `align-and-count.sh` per study, `04_bulk_te_deseq.R`
   per study, `04 baseline` and `04 compare`.

Before trusting any study, check the perturbed gene's own log2FC in its DE
table. Knockdowns that did not knock down are common.

## 6. Caveats

- **Context mismatch.** Most bulk perturbations are cell lines (K562,
  HEK293, HepG2, mESCs) with complete or acute loss. ISP runs in primary
  single cells from another tissue. Disagreement is evidence about the
  model *in that context*, not a verdict on the biology.
- **Multi-mapping.** TEcount `--mode multi` distributes multimappers by EM.
  Young families (L1HS, AluY, SVA) carry the most uncertainty. Unique-read
  sensitivity, as done in TE_Ocean_Acidification
  (`te-uniqueread-sensitivity-20260903`), is the natural robustness check.
- **Read length and library type** vary across studies. Strandedness is
  taken from metadata when given. The inference PerturbAtlas uses (STAR
  ReadsPerGene, 80% rule) can be reused when it is not.
- **Selection bias.** Published perturbation studies over-represent
  regulators already known to matter, so they favour positive controls and
  say little about novel ISP hits.
- **Fish.** PerturbAtlas includes zebrafish. That opens a route to checking
  teleost TE responses at the regulator level, which is relevant to the
  *A. polyacanthus* work. No family-level transfer is possible.
