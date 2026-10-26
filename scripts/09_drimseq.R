#!/usr/bin/env Rscript
# Test differential poly(A) site usage between stages with DRIMSeq.
#
# Input: the folder written by step 08 (site_counts.tsv, samples.tsv).
# Output in the same folder:
#   drimseq_gene.tsv      one row per tested gene: DRIMSeq p value and BH adjusted
#                         p value, the stageR gene adjusted p value (gene_padj),
#                         gene_confirmed, and fit_note
#   drimseq_site.tsv      one row per site: DRIMSeq site p value, fitted usage per
#                         stage, the stageR site adjusted p value (site_padj) and
#                         site_confirmed
#   drimseq_stager.tsv    the stageR table alone
#   untestable_genes.tsv  genes left out and why (depth, one observed site, fit failed)
#
# The stageR decisions (gene_confirmed, site_confirmed) are the ones every later
# step uses. They control the overall false discovery rate across the gene
# screen and the site confirmation together.
#
# Usage: Rscript 09_drimseq.R <folder> [fdr] [min_reads_per_gene] [min_reads_per_site] [pseudocount] [fallback_pseudocount]
#
# Counts are fitted as observed (pseudocount 0). Features with no read in any
# sample are dropped, and a gene with a single observed feature is untestable
# for relative usage. A site with zero reads in every replicate of one stage can
# make the Dirichlet-multinomial fit fail, which DRIMSeq reports as a missing p
# value. Those genes alone are refitted with the fallback pseudocount (default
# 0.5) added to every cell of that gene, and fit_note records it.
suppressPackageStartupMessages({library(DRIMSeq); library(stageR)})
args <- commandArgs(trailingOnly = TRUE)
folder <- args[1]
fdr <- if (length(args) >= 2) as.numeric(args[2]) else 0.05
min_gene <- if (length(args) >= 3) as.integer(args[3]) else 10
min_site <- if (length(args) >= 4) as.integer(args[4]) else 0
pseudo <- if (length(args) >= 5) as.numeric(args[5]) else 0
fallback <- if (length(args) >= 6) as.numeric(args[6]) else 0.5

counts <- read.delim(file.path(folder, "site_counts.tsv"), check.names = FALSE)
# DRIMSeq needs feature names that are unique across genes
counts$feature_id <- paste(counts$gene_id, counts$feature_id, sep = "|")
samples <- read.delim(file.path(folder, "samples.tsv"))
samples$stage <- factor(samples$stage, levels = c("BSF", "PCF"))
count_cols <- samples$sample_id
untestable <- data.frame(gene_id = character(), reason = character())

# Depth filter on the raw counts: genes with at least min_gene reads in every
# replicate; sites with at least min_site reads in three replicates.
gene_depth <- aggregate(counts[count_cols], by = list(gene_id = counts$gene_id), FUN = sum)
deep <- gene_depth$gene_id[apply(gene_depth[count_cols], 1, min) >= min_gene]
shallow <- setdiff(gene_depth$gene_id, deep)
if (length(shallow)) untestable <- rbind(untestable, data.frame(gene_id = shallow,
  reason = sprintf("fewer than %d reads in at least one replicate", min_gene)))
counts <- counts[counts$gene_id %in% deep & rowSums(counts[count_cols] >= min_site) >= 3, ]
# Features with no read anywhere carry no information and distort proportions once a pseudocount is added
counts <- counts[rowSums(counts[count_cols]) > 0, ]
n_feat <- table(counts$gene_id)
single <- names(n_feat)[n_feat < 2]
if (length(single)) untestable <- rbind(untestable, data.frame(gene_id = single,
  reason = "one observed site only: relative usage cannot change"))
counts <- counts[!counts$gene_id %in% single, ]

fit_gene_set <- function(cnt, add) {
  cnt[count_cols] <- cnt[count_cols] + add
  d <- dmDSdata(counts = cnt, samples = samples)
  d <- dmFilter(d)
  design <- model.matrix(~ stage, data = samples)
  set.seed(42)
  d <- dmPrecision(d, design = design)
  d <- dmFit(d, design = design)
  d <- dmTest(d, coef = "stagePCF")
  gene_res <- results(d)
  site_res <- results(d, level = "feature")
  prop <- proportions(d)
  bsf <- rowMeans(prop[, samples$sample_id[samples$stage == "BSF"], drop = FALSE])
  pcf <- rowMeans(prop[, samples$sample_id[samples$stage == "PCF"], drop = FALSE])
  site_res$usage_BSF <- round(bsf[match(site_res$feature_id, prop$feature_id)], 3)
  site_res$usage_PCF <- round(pcf[match(site_res$feature_id, prop$feature_id)], 3)
  list(gene = gene_res, site = site_res)
}

main <- fit_gene_set(counts, pseudo)
main$gene$fit_note <- if (pseudo == 0) "observed counts" else sprintf("pseudocount %g", pseudo)
failed <- main$gene$gene_id[is.na(main$gene$pvalue)]
cat(sprintf("%d genes fitted on observed counts; %d without a p value\n", nrow(main$gene), length(failed)))
if (length(failed) && fallback > 0) {
  fb <- fit_gene_set(counts[counts$gene_id %in% failed, ], fallback)
  fb$gene$fit_note <- sprintf("pseudocount %g added after the fit on observed counts failed (zero counts in one stage)", fallback)
  still <- fb$gene$gene_id[is.na(fb$gene$pvalue)]
  main$gene <- rbind(main$gene[!main$gene$gene_id %in% failed, ], fb$gene)
  main$site <- rbind(main$site[!main$site$gene_id %in% failed, ], fb$site)
  cat(sprintf("%d of them refitted with pseudocount %g; %d still without a p value\n",
              length(failed), fallback, length(still)))
  failed <- still
}
if (length(failed)) untestable <- rbind(untestable, data.frame(gene_id = failed, reason = "model fit failed"))
gene_res <- main$gene[!is.na(main$gene$pvalue), ]
site_res <- main$site[main$site$gene_id %in% gene_res$gene_id, ]
# BH over the genes that were actually tested
gene_res$adj_pvalue <- p.adjust(gene_res$pvalue, method = "BH")
site_res$usage_change_PCF_minus_BSF <- site_res$usage_PCF - site_res$usage_BSF

# stageR: screen on the gene p values, confirm sites, overall FDR = fdr.
# Every feature of a tested gene stays in the confirmation family. A feature
# without a p value (fit failed for that feature) is kept with p = 1, a
# non-discovery, so that the two-feature rule of the DTU adjustment is not
# applied to a composition the gene test did not see. Nothing is removed after
# fitting.
site_ok <- site_res
site_ok$pvalue[is.na(site_ok$pvalue)] <- 1
site_res$site_p_missing <- is.na(site_res$pvalue)
n_sites <- table(site_ok$gene_id)
keep <- names(n_sites)[n_sites >= 2]
site_ok <- site_ok[site_ok$gene_id %in% keep, ]
gene_ok <- gene_res[gene_res$gene_id %in% keep, ]
# stageR splits identifiers on ':', so it gets cleaned keys and the result is mapped back
safe <- function(x) gsub("[: ]", "_", x)
back_gene <- setNames(gene_ok$gene_id, safe(gene_ok$gene_id)); back_feat <- setNames(site_ok$feature_id, safe(site_ok$feature_id))
pScreen <- gene_ok$pvalue; names(pScreen) <- safe(gene_ok$gene_id)
pConf <- matrix(site_ok$pvalue, ncol = 1, dimnames = list(safe(site_ok$feature_id), "stagePCF"))
tx2gene <- data.frame(feature_id = safe(site_ok$feature_id), gene_id = safe(site_ok$gene_id))
st <- stageRTx(pScreen = pScreen, pConfirmation = pConf, pScreenAdjusted = FALSE, tx2gene = tx2gene)
st <- stageWiseAdjustment(st, method = "dtu", alpha = fdr, allowNA = TRUE)
padj <- getAdjustedPValues(st, order = TRUE, onlySignificantGenes = FALSE)
colnames(padj) <- c("gene_id", "feature_id", "gene_padj", "site_padj")
padj$gene_id <- unname(back_gene[padj$gene_id]); padj$feature_id <- unname(back_feat[padj$feature_id])
padj$gene_confirmed <- padj$gene_padj <= fdr
padj$site_confirmed <- !is.na(padj$site_padj) & padj$site_padj <= fdr

gene_res$gene_padj <- padj$gene_padj[match(gene_res$gene_id, padj$gene_id)]
gene_res$gene_confirmed <- gene_res$gene_id %in% padj$gene_id[padj$gene_confirmed]
site_res$site_padj <- padj$site_padj[match(site_res$feature_id, padj$feature_id)]
site_res$site_confirmed <- site_res$feature_id %in% padj$feature_id[padj$site_confirmed]
strip <- function(x) sub("^[^|]*\\|", "", x)
site_res$feature_id <- strip(site_res$feature_id)
padj$feature_id <- strip(padj$feature_id)

write.table(gene_res[order(gene_res$pvalue), ], file.path(folder, "drimseq_gene.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(site_res[order(site_res$pvalue), ], file.path(folder, "drimseq_site.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(padj, file.path(folder, "drimseq_stager.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(untestable, file.path(folder, "untestable_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("stageR: %d of %d tested genes confirmed at overall FDR %.2f, %d sites confirmed; %d genes untestable; %d site p values missing, kept as 1\n",
            sum(gene_res$gene_confirmed), nrow(gene_res), fdr, sum(site_res$site_confirmed), nrow(untestable), sum(site_res$site_p_missing)))
