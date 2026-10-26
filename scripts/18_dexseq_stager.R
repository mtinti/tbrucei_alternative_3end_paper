#!/usr/bin/env Rscript
# Test differential end-cluster usage between stages with DEXSeq, then stageR (docs/methods.md, section 3).
#
# Input: a folder from step 17 with site_counts.tsv (gene_id, feature_id, one
# column per replicate; supported clusters and "other" of the testable windows)
# and samples.tsv (sample_id, stage).
#
# Model, as DEXSeq defines it: for every feature the counts of this feature
# against the rest of its window; full design ~ sample + exon + stage:exon,
# reduced design ~ sample + exon, where exon is DEXSeq's this-feature-against-
# the-rest factor. Dispersions are per feature, shrunk toward the trend.
# "other" is a feature like the clusters. DEXSeq's independent filtering is off
# because steps 16 and 17 already filtered. A feature whose test gives no p
# value (fit failure) is kept with p = 1 and flagged; nothing is removed after
# fitting.
#
# Gene level: DEXSeq's per-gene q value screens the genes. It is already
# adjusted, so stageR takes it with pScreenAdjusted = TRUE; the features'
# unadjusted p values are the confirmation stage. Overall false discovery
# rate = fdr, gene-level OFDR.
#
# Output in the folder: dexseq_site.tsv (feature p value, p_missing, DEXSeq
# log2 fold change of relative usage, raw usage fraction per replicate, mean
# raw fraction per stage and their difference, stageR site_padj and
# site_confirmed), dexseq_gene.tsv (per-gene q value, gene_confirmed,
# features, failed_features), stager.tsv.
#
# Usage: Rscript 18_dexseq_stager.R <folder> [fdr] [threads]
suppressPackageStartupMessages({library(DEXSeq); library(stageR); library(BiocParallel)})
args <- commandArgs(trailingOnly = TRUE)
folder <- args[1]
fdr <- if (length(args) >= 2) as.numeric(args[2]) else 0.05
threads <- if (length(args) >= 3) as.integer(args[3]) else 4
BPPARAM <- MulticoreParam(threads)

counts <- read.delim(file.path(folder, "site_counts.tsv"), check.names = FALSE)
samples <- read.delim(file.path(folder, "samples.tsv"))
samples$stage <- factor(samples$stage, levels = c("BSF", "PCF"))
reps <- samples$sample_id
feature_key <- paste(counts$gene_id, counts$feature_id, sep = "|")
# DEXSeq removes ':' and spaces from identifiers; give it clean keys and map back afterwards
safe <- function(x) gsub("[: ]", "_", x)
safe_gene <- safe(counts$gene_id); safe_key <- safe(feature_key)
stopifnot(!anyDuplicated(safe_key))
back_gene <- setNames(counts$gene_id, safe_gene); back_key <- setNames(feature_key, safe_key)
mat <- as.matrix(counts[reps]); rownames(mat) <- feature_key
sample_data <- data.frame(row.names = reps, stage = samples$stage)

dxd <- DEXSeqDataSet(countData = mat, sampleData = sample_data,
                     design = ~ sample + exon + stage:exon,
                     featureID = safe_key, groupID = safe_gene)
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd, BPPARAM = BPPARAM)
dxd <- testForDEU(dxd, reducedModel = ~ sample + exon, fullModel = ~ sample + exon + stage:exon, BPPARAM = BPPARAM)
dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "stage", BPPARAM = BPPARAM)
res <- DEXSeqResults(dxd, independentFiltering = FALSE)
res_df <- as.data.frame(res[, c("groupID", "featureID", "pvalue", "log2fold_PCF_BSF")])
colnames(res_df) <- c("gene_id", "feature_key", "pvalue", "log2fold_usage_PCF_vs_BSF")
res_df$gene_id <- unname(back_gene[as.character(res_df$gene_id)])
res_df$feature_key <- unname(back_key[as.character(res_df$feature_key)])
stopifnot(!anyNA(res_df$gene_id), !anyNA(res_df$feature_key))
res_df$feature_id <- sub("^[^|]*\\|", "", res_df$feature_key)
res_df$p_missing <- is.na(res_df$pvalue)
cat(sprintf("%d features in %d genes; %d features without a p value, kept with p = 1\n",
            nrow(res_df), length(unique(res_df$gene_id)), sum(res_df$p_missing)))
res_df$pvalue[res_df$p_missing] <- 1

# raw usage fractions per replicate and per stage, the effect size that is reported
frac <- sweep(mat, 2, 1, "/")
for (g in unique(counts$gene_id)) {
  idx <- counts$gene_id == g
  tot <- colSums(mat[idx, , drop = FALSE])
  frac[idx, ] <- sweep(mat[idx, , drop = FALSE], 2, pmax(tot, 1), "/")
}
frac <- frac[res_df$feature_key, , drop = FALSE]
colnames(frac) <- paste0("frac_", reps)
res_df <- cbind(res_df, round(frac, 4))
res_df$usage_BSF <- round(rowMeans(frac[, paste0("frac_", reps[samples$stage == "BSF"]), drop = FALSE]), 4)
res_df$usage_PCF <- round(rowMeans(frac[, paste0("frac_", reps[samples$stage == "PCF"]), drop = FALSE]), 4)
res_df$usage_change_points <- round(100 * (res_df$usage_PCF - res_df$usage_BSF), 1)

# gene level: DEXSeq per-gene q value over the features with their (p = 1 substituted) p values
res_for_q <- res
res_for_q@listData[["pvalue"]] <- res_df$pvalue[match(as.character(res$featureID), safe(res_df$feature_key))]
stopifnot(!anyNA(res_for_q@listData[["pvalue"]]))
qgene <- perGeneQValue(res_for_q)
gene_df <- data.frame(gene_id = unname(back_gene[names(qgene)]), gene_qvalue = as.numeric(qgene))
gene_df$features <- as.integer(table(res_df$gene_id)[gene_df$gene_id])
gene_df$failed_features <- as.integer(tapply(res_df$p_missing, res_df$gene_id, sum)[gene_df$gene_id])
missing_genes <- setdiff(unique(res_df$gene_id), gene_df$gene_id)
if (length(missing_genes)) {
  gene_df <- rbind(gene_df, data.frame(gene_id = missing_genes, gene_qvalue = NA,
                                       features = as.integer(table(res_df$gene_id)[missing_genes]),
                                       failed_features = as.integer(tapply(res_df$p_missing, res_df$gene_id, sum)[missing_genes])))
}
gene_df$screen_outcome <- ifelse(is.na(gene_df$gene_qvalue), "no gene q value (all features failed)", "screened")

# stageR: screen on the adjusted gene q values, confirm features on their raw p values
# stageR splits identifiers on ':' as well, so it gets the clean keys and the result is mapped back
ok <- gene_df[!is.na(gene_df$gene_qvalue), ]
pScreen <- ok$gene_qvalue; names(pScreen) <- safe(ok$gene_id)
conf <- res_df[res_df$gene_id %in% ok$gene_id, ]
pConf <- matrix(conf$pvalue, ncol = 1, dimnames = list(safe(conf$feature_key), "stage"))
tx2gene <- data.frame(feature_key = safe(conf$feature_key), gene_id = safe(conf$gene_id))
st <- stageRTx(pScreen = pScreen, pConfirmation = pConf, pScreenAdjusted = TRUE, tx2gene = tx2gene)
st <- stageWiseAdjustment(st, method = "dtu", alpha = fdr, allowNA = TRUE)
padj <- getAdjustedPValues(st, order = TRUE, onlySignificantGenes = FALSE)
colnames(padj) <- c("gene_id", "feature_key", "gene_padj", "site_padj")
padj$gene_id <- unname(back_gene[padj$gene_id]); padj$feature_key <- unname(back_key[padj$feature_key])
stopifnot(!anyNA(padj$gene_id), !anyNA(padj$feature_key))
padj$gene_confirmed <- padj$gene_padj <= fdr
padj$site_confirmed <- !is.na(padj$site_padj) & padj$site_padj <= fdr

gene_df$gene_padj <- padj$gene_padj[match(gene_df$gene_id, padj$gene_id)]
gene_df$gene_confirmed <- gene_df$gene_id %in% padj$gene_id[padj$gene_confirmed]
res_df$site_padj <- padj$site_padj[match(res_df$feature_key, padj$feature_key)]
res_df$site_confirmed <- res_df$feature_key %in% padj$feature_key[padj$site_confirmed]
padj$feature_id <- sub("^[^|]*\\|", "", padj$feature_key)

out_cols <- c("gene_id", "feature_id", "pvalue", "p_missing", "log2fold_usage_PCF_vs_BSF", paste0("frac_", reps),
              "usage_BSF", "usage_PCF", "usage_change_points", "site_padj", "site_confirmed")
write.table(res_df[order(res_df$pvalue), out_cols], file.path(folder, "dexseq_site.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gene_df[order(gene_df$gene_qvalue), ], file.path(folder, "dexseq_gene.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(padj[, c("gene_id", "feature_id", "gene_padj", "site_padj", "gene_confirmed", "site_confirmed")],
            file.path(folder, "stager.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("stageR: %d of %d screened genes confirmed at gene-level OFDR %.2f; %d features confirmed; %d genes without a gene q value\n",
            sum(gene_df$gene_confirmed), sum(!is.na(gene_df$gene_qvalue)), fdr, sum(res_df$site_confirmed),
            sum(is.na(gene_df$gene_qvalue))))
