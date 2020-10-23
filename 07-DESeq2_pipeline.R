# This script was used to call differentially expressed genes in cold-treated samples on both plaNET-Seq and TSS-Seq data;
# It also contains code to reproduce figures 4c and S5a-b;

library(DESeq2)
library(rtracklayer)
library(ggplot2)
library(reshape2)

r_dir <- "." # change to the folder where you saved the custom functions from https://github.com/Maxim-Ivanov/Kindgren_et_al_2019
scripts <- c("batchReadTrackData.R", "getOverlappingScores.R")
for (script in scripts) { source(file.path(r_dir, script)) }

genes_araport_adj <- readRDS("genes_araport_adj.RDS")
novel <- readRDS("Novel_groHMM_transcripts.RDS")

txdb_araport <- makeTxDbFromGFF("Araport11_GFF3_genes_transposons.201606.gff.gz")
ebg_araport <- exonsBy(txdb_araport, by = "gene") # exons grouped by gene
seqinfo(ebg_araport, new2old = as.integer(c(1:5, 7, 6))) <- seqinfo(genes_araport_adj)

### Call differential transcription on plaNET-Seq data:

planet_dir <- "." # change to the directory containing plaNET-Seq Bedgraph files obtained from 02-Postprocessing_plaNET-Seq.R
planet_files <- list.files(planet_dir, pattern = "biorep.*fw_rev.bedgraph.gz$")
# (check the number and order of elements in <planet_files>. Ensure that it fits the <coldata> dataframe created below. Subset and/or reorder <planet_files> if necessary)
planet_data <- batchReadTrackData(planet_files, dir = planet_dir, format = "bedGraph", seqinfo = seqinfo(genes_araport_adj))
names(planet_data) <- paste0("plaNET_", sub("_fw_rev.bedgraph.gz", "", names(planet_data)))

genes_cm <- getOverlappingScores(genes_araport_adj, planet_data, value = "count_matrix")
novel_cm <- getOverlappingScores(novel, planet_data, value = "count_matrix")
coldata <- data.frame("Sample" = rep(c("WT", "Cold_3h", "Cold_12h", "Fas2", "DMSO", "PlaB"), each = 2),
                      "Genotype" = rep(c("wt", "fas2", "wt"), c(6, 2, 4)),
                      "Treatment" = rep(c("No", "Cold_3h", "Cold_12h", "No", "DMSO", "PlaB"), each = 2),
                      "Replicate" = paste0("Rep", 1:2),
                      row.names = names(planet_data))
genes_se <- SummarizedExperiment(rowRanges = genes_araport_adj, colData = coldata, assays = list(counts = genes_cm))
novel_se <- SummarizedExperiment(rowRanges = novel, colData = coldata, assays = list(counts = novel_cm))

deseq_pipeline <- function(se, lfc = 1, pval = 0.05) {
  dds <- DESeqDataSet(se, design = ~ Sample)
  dds$Sample <- relevel(dds$Sample, ref = "WT")
  dds <- DESeq(dds)
  contrasts <- list(c("Cold_3h", "WT"), c("Cold_12h", "WT"), c("Cold_12h", "Cold_3h"), c("Fas2", "WT"), c("PlaB", "DMSO"))
  all_res <- vector("list", length(contrasts))
  for (i in seq_along(contrasts)) {
    first <- contrasts[[i]][[1]]
    second <- contrasts[[i]][[2]]
    res <- results(dds, contrast = c("Sample", first, second))
    res <- lfcShrink(dds, contrast = c("Sample", first, second), res = res, type = "normal")
    res_out <- res[c("log2FoldChange", "padj")]
    res_out$de <- "No"
    res_out$de[res_out$log2FoldChange >= lfc & res_out$padj <= pval] <- "Up"
    res_out$de[res_out$log2FoldChange <= -lfc & res_out$padj <= pval] <- "Down"
    all_res[[i]] <- res_out
  }
  df <- as.data.frame(do.call(cbind, all_res))
  samples <- unlist(lapply(contrasts, paste, collapse = "_vs_"))
  lfc <- paste0(samples, "_lfc")
  padj <- paste0(samples, "_padj")
  de <- paste0(samples, "_de")
  names(df) <- as.character(mapply(function(x, y, z) {c(x, y, z)}, lfc, padj, de))
  mcols(se) <- cbind(mcols(se), df)
  return(se)
}

genes_de <- deseq_pipeline(genes_se)
saveRDS(genes_de, "genes_de.RDS")
novel_de <- deseq_pipeline(novel_se)
saveRDS(novel_de, "novel_de.RDS")

### Fig. S5B (number of exons in DE genes):

# Extract exon numbering within each gene:
ebg <- exonsBy(txdb, by = "gene")
exon_names <- mcols(unlist(ebg))$exon_name
exon_nums <- unlist(lapply(strsplit(exon_names, ".", fixed = TRUE), function(x) { sub("exon", "", x[[3]]) }))
exons_unl <- unlist(ebg)
mcols(exons_unl)$Num <- as.integer(exon_nums)
ebg <- relist(exons_unl, ebg)

genes_gr <- rowRanges(genes_de)
genes_pcd <- genes_gr[mcols(genes_gr)$tx_type == "mRNA" & seqnames(genes_gr) %in% 1:5]
all_cols <- names(mcols(genes_pcd))
de_cols <- all_cols[grepl("_de", all_cols)]

# Print exon count distribution:
res <- vector("list", length(de_cols) * 3)
for (i in seq_along(de_cols)) {
  de_col <- de_cols[[i]]
  decisions <- mcols(genes3_pcd)[, de_col]
  for (j in 1:3) {
    tag <- c("Up", "Down", "No")[[j]]
    gene_names <- names(genes3_pcd)[decisions == tag]
    exon_counts <- lengths(ebg[names(ebg) %in% gene_names])
    report <- vector("integer", 11)
    for (k in 1:10) {
      report[[k]] <- sum(exon_counts == k)
    }
    report[[11]] <- sum(exon_counts > 10)
    res[[(i - 1) * 3 + j]] <- c(list(de_col, tag), report)
  }
}
res_df <- as.data.frame(do.call(rbind, res))
for (i in 1:ncol(res_df)) { res_df[[i]] <- unlist(res_df[[i]]) }
names(res_df) <- c("Contrast", "DE", 1:11)
mat <- as.matrix(res_df[, 3:13])
res_df[, 3:13] <- mat / rowSums(mat)

long_df <- melt(res_df, id.vars = 1:2, variable.name = "Exon_count", value.name = "N_genes")
long_df$Exon_count <- as.integer(long_df$Exon_count)
long_df$Contrast <- factor(long_df$Contrast)
long_df$DE <- factor(long_df$DE)

ttl <- "Exon counts in DE genes (Fig. S5b)"
p <- ggplot(long_df, aes(x = Exon_count, y = N_genes, fill = DE)) +
  geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.5) +
  facet_grid(Contrast ~ .) + ggtitle(ttl)
for (ext in c(".png", ".pdf")) {
  ggsave(paste0(ttl, ext), plot = p, width = 6, height = 10, units = "in")
}


### Fig. S5a (length of Cold DE genes):

out <- vector("list", length(de_cols))
for (i in seq_along(de_cols)) {
  de_col <- de_cols[[i]]
  df <- data.frame("Length" = width(genes_pcd), "DE" = mcols(genes_pcd)[, de_col], "Group" = sub("_de", "", de_col))
  out[[i]] <- df
}
results <- do.call(rbind, out)

p <- ggplot(results, aes(x = DE, y = Length)) + geom_boxplot(outlier.colour = NA) + ylim(0, 6000) + facet_grid(. ~ Group)
for (ext in c(".png", ".pdf")) {
  ggsave(paste0("Length of DE genes (Fig. S5a)", ext), plot = p, width = 12, height = 8, units = "in")
}

### Call differential expression on TSS-Seq data:

# Load TSS-Seq BigWig files:
tss_dir <- "." # change to the directory where you downloaded all BigWig files from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119304 (Kindgren at al., 2018 - PMID 30385760)
tss_data_fw <- batchReadTrackData(list.files(tss_dir, pattern = "Plus.bw$"), dir = tss_dir, format = "BigWig")
strand(tss_data_fw) <- "+"
tss_data_rev <- batchReadTrackData(list.files(tss_dir, pattern = "Minus.bw$"), dir = tss_dir, format = "BigWig")
strand(tss_data_rev) <- "-"
# Merge Fw and Rev files:
tss_data <- mapply(function(x, y) { return(sort(c(x, y))) }, tss_data_fw, tss_rev, SIMPLIFY = FALSE)
names(tss_data) <- c("RT_rep1", "RT_rep2", "Cold_rep1", "Cold_rep2")

# Quantify TSS-Seq data on gene promoters:
genes_prom <- suppressWarnings(trim(resize(resize(genes_araport_adj, 300, "start"), 500, "end"))) # [TSS-200, TSS+300]
genes_tss_cm <- getOverlappingScores(genes_prom, tss_data, value = "count_matrix")

# Make SummarizedExperiment objects:
coldata <- data.frame("Sample" = paste("s", 1:4),
                      "Genotype" = "wt",
                      "Treatment" = rep(c("RT", "Cold_3h"), each = 2),
                      "Replicate" = paste0("rep", 1:2),
                      row.names = colnames(genes_tss_cm))
genes_tss_se <- SummarizedExperiment(rowRanges = genes_prom, colData = coldata, assays = list(counts = genes_tss_cm))

tss_deseq_pipeline <- function(se) {
  dds <- DESeqDataSet(se, design = ~ Sample)
  dds$Treatment <- relevel(dds$Treatment, ref = "RT")
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("Treatment", "Cold_3h", "RT"))
  res <- lfcShrink(dds, contrast = c("Treatment", "Cold_3h", "RT"), res = res)
  de <- rep("No", nrow(res))
  de[res$log2FoldChange >= 1 & res2$padj <= 0.05] <- "Up"
  de[res$log2FoldChange <= -1 & res2$padj <= 0.05] <- "Down"
  return(de)
}

genes_tss_decisions <- tss_deseq_pipeline(genes_tss_se)

# Fig. 4c (comparison of DE decisions for the same genes on plaNET-Seq vs TSS-Seq data):
table(mcols(genes_de)$Cold_3h_vs_WT_de, genes3_tss_decisions)


