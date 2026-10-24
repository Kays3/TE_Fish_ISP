#!/usr/bin/env Rscript
## Per-study TE differential expression (perturbed vs control) from TEcount
## .cntTable files, for 04_bulk_validation.py. Mirrors the DESeq2 setup in
## TE_Ocean_Acidification/pipelines/shared/deseq-te-counts.R.
##
## Usage: Rscript 04_bulk_te_deseq.R COUNTS_DIR SAMPLESHEET_TSV OUTPUT_TSV

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript 04_bulk_te_deseq.R COUNTS_DIR SAMPLESHEET_TSV OUTPUT_TSV")
}
suppressPackageStartupMessages({
  library(DESeq2)
})

counts_dir <- normalizePath(args[[1]], mustWork = TRUE)
sheet <- read.delim(args[[2]], stringsAsFactors = FALSE)
out_path <- args[[3]]

files <- list.files(counts_dir, pattern = "\\.cntTable$", full.names = TRUE, recursive = TRUE)
if (!length(files)) stop("No TEcount .cntTable files under ", counts_dir)

read_one <- function(path) {
  tab <- read.delim(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  setNames(as.numeric(tab[[2]]), tab[[1]])
}
tabs <- lapply(files, read_one)
names(tabs) <- sub("\\.cntTable$", "", basename(files))
features <- unique(unlist(lapply(tabs, names)))
mat <- sapply(tabs, function(v) { x <- v[features]; x[is.na(x)] <- 0; x })
rownames(mat) <- features

keep <- intersect(sheet$sample_id, colnames(mat))
if (length(keep) != nrow(sheet)) {
  stop("Missing counts for: ", paste(setdiff(sheet$sample_id, colnames(mat)), collapse = ", "))
}
sheet <- sheet[match(keep, sheet$sample_id), ]
coldata <- data.frame(role = factor(sheet$treatment, levels = c("control", "perturbed")),
                      row.names = sheet$sample_id)
## Optional pairing/blocking factor (e.g. clutch): design ~ block + role.
## Used only when it has >1 level and every level appears in both roles or
## at least one block spans both roles (otherwise it is confounded with role).
design <- ~ role
if ("block" %in% names(sheet) && any(nzchar(sheet$block))) {
  coldata$block <- factor(sheet$block)
  spans <- tapply(coldata$role, coldata$block, function(x) length(unique(x)) == 2)
  if (nlevels(coldata$block) > 1 && any(spans)) design <- ~ block + role
}
message("design: ", deparse(design))

dds <- DESeqDataSetFromMatrix(round(mat[, keep]), coldata, design = design)
min_reps <- min(table(coldata$role))
dds <- dds[rowSums(counts(dds) >= 10L) >= min_reps, ]
dds <- DESeq(dds)
res <- results(dds, contrast = c("role", "perturbed", "control"))
res <- lfcShrink(dds, coef = "role_perturbed_vs_control", res = res, type = "apeglm")

out <- data.frame(feature = rownames(res), baseMean = res$baseMean,
                  log2FoldChange = res$log2FoldChange, padj = res$padj)
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write.table(out, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), sub("\\.tsv$", ".sessionInfo.txt", out_path))
