#!/usr/bin/env Rscript
# Test the abundance of the independent transcripts between stages with edgeR (docs/methods.md, section 5).
#
# Input: the step 22 folder (catalogue.tsv with full-length read counts per
# replicate, gene_full_length_counts.tsv) and the samples table.
# Normalisation: TMM effective library sizes computed from the full-length
# read counts of all protein-coding genes and all catalogue records together,
# never from the hits. Test: quasi-likelihood F test, stage the only factor,
# tagwise dispersion, BH over the catalogue records that pass the expression
# filter. The result is a change in relative abundance within the full-length
# read population.
# Output: transcript_abundance.tsv in the folder.
# Usage: Rscript 23_transcript_abundance.R <folder> <samples.tsv>
suppressPackageStartupMessages(library(edgeR))
args <- commandArgs(trailingOnly = TRUE)
folder <- args[1]; samples <- read.delim(args[2])
samples$stage <- factor(samples$stage, levels = c("BSF", "PCF")); reps <- samples$sample_id
cat_tab <- read.delim(file.path(folder, "catalogue.tsv"), check.names = FALSE)
genes <- read.delim(file.path(folder, "gene_full_length_counts.tsv"), check.names = FALSE)
ref <- rbind(as.matrix(genes[reps]), as.matrix(cat_tab[reps]))
rownames(ref) <- c(paste0("gene|", genes$gene_id), paste0("record|", cat_tab$record))
y <- DGEList(counts = ref, group = samples$stage)
y <- calcNormFactors(y, method = "TMM")
cat(sprintf("TMM factors from %d genes and %d records: %s\n", nrow(genes), nrow(cat_tab), paste(round(y$samples$norm.factors, 3), collapse = " ")))
keep <- rownames(ref) %in% paste0("record|", cat_tab$record[cat_tab$in_family == "True" | cat_tab$in_family == TRUE])
yt <- y[keep, , keep.lib.sizes = TRUE]     # library sizes and norm factors stay those of the whole reference
design <- model.matrix(~ stage, data = samples)
yt <- estimateDisp(yt, design)
fit <- glmQLFit(yt, design)
res <- glmQLFTest(fit, coef = "stagePCF")
tab <- topTags(res, n = Inf, sort.by = "none")$table
tab$record <- sub("^record\\|", "", rownames(tab))
cpm_mat <- cpm(y)[rownames(tab), , drop = FALSE]
tab$mean_cpm_BSF <- round(rowMeans(cpm_mat[, reps[samples$stage == "BSF"], drop = FALSE]), 1)
tab$mean_cpm_PCF <- round(rowMeans(cpm_mat[, reps[samples$stage == "PCF"], drop = FALSE]), 1)
out <- merge(cat_tab[c("record", "kind", "members", "chrom", "start", "end", "strand", "annotated", reps)],
             tab[c("record", "logFC", "logCPM", "F", "PValue", "FDR", "mean_cpm_BSF", "mean_cpm_PCF")], by = "record")
colnames(out)[colnames(out) == "logFC"] <- "log2FC_PCF_vs_BSF"
out$abundance_changed <- out$FDR <= 0.05
write.table(out[order(out$PValue), ], file.path(folder, "transcript_abundance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("abundance: %d of %d records tested change at BH FDR 0.05\n", sum(out$abundance_changed), nrow(out)))
