# This code was used to filter BAM files produced by https://github.com/Maxim-Ivanov/Kindgren_et_al_2019/01-Alignment_plaNET-Seq.sh
# It returns Bedgraph files which are used by the downstream scripts in this study;
# It also returns strand-specific BigWig files which were uploaded to GEO;
# Moreover, it allows to reproduce Fig. S1D;

library(GenomicAlignments)
library(dplyr)
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28
ebg1 <- exonsBy(txdb, by = "gene")

txdb3 <- makeTxDbFromGFF("Araport11_GFF3_genes_transposons.201606.gff.gz")
ebg3 <- exonsBy(txdb3, by = "gene")
seqinfo(ebg3, new2old = as.integer(c(1:5, 7, 6))) <- seqinfo(txdb)

library(BSgenome.Athaliana.TAIR.TAIR9)
bs <- BSgenome.Athaliana.TAIR.TAIR9
seqlevels(bs) <- seqlevels(txdb)

scripts <- c("convertGAlignmentsToCoverage.R", "saveGRangesAsBedGraph.R", "mergeBedgraphs.R", "loadNETSeqBAM.R", 
             "removeFirstAndLastExons.R", "getOverlappingScores.R", "normalizeGR.R")
for (script in scripts) { source(script) }

findSpliceSites <- function(ebg) {
  res <- removeFirstAndLastExons(ebg)
  exons_wo_first <- res[[1]]
  exons_wo_last <- res[[2]]
  ss_5p <- GenomicRanges::reduce(resize(resize(exons_wo_last, 1, "end"), 3, "center")) # the last base of exon (+ 1 base up  + 1 base down)
  ss_3p <- GenomicRanges::reduce(resize(flank(exons_wo_first, 1), 3, "center")) # the last base of the upstream intron (+ 1 base up  + 1 base down)
  return(list(ss_5p, ss_3p))
}

ss1 <- findSpliceSites(ebg1)
donor_1 <- ss1[[1]]
acc_1 <- ss1[[2]]
ss2 <- findSpliceSites(ebg3)
donor_2 <- ss2[[1]]
acc_2 <- ss2[[2]]
all_donor <- sort(GenomicRanges::reduce(c(donor_1, donor_2))) # all donor splice sites
all_acc <- sort(GenomicRanges::reduce(c(acc_1, acc_2))) # all acceptor splice sites

bamdir <- "." # change to path of your folder which contains the BAM files returned by 01-Alignment_plaNET-Seq.sh
bamfiles <- list.files(bamdir, pattern = "full.*mapq.bam$")

# Load and filter all BAM files, save sequencing coverage as Bedgraph files:
# (the Bedgraph files contain both positive and negative values in the 4th column to encode the strand information)
df <- loadNETSeqBAM(bamfiles, bamdir, all_donor, all_acc, bs, trim_names = "_full_sorted_dedup_clean_mapq")

# Merge biological replicates (without normalization to 1M reads):
rep1 <- list.files(bamdir, pattern = "biorep1.*bedgraph.gz$")
rep2 <- sub("biorep1", "biorep2", rep1)

for (i in seq_along(rep1)) {
  f1 <- rep1[[i]]
  f2 <- rep2[[i]]
  message(f1, " + ", f2)
  gr <- mergeBedgraphs(bg_fw_rev = file.path(bamdir, list(f1, f2)), seqinfo = seqinfo(txdb), norm = FALSE)
  saveGRangesAsBedGraph(gr, file.path(bamdir, sub("biorep1", "merged", f1)))
}

# Convert Bedgraph to BigWig (for GEO):
bg_files <- list.files(bamdir, pattern = "bedgraph.gz$")
for (i in seq_along(bg_files)) {
  bg <- bg_files[[i]]
  message(bg)
  gr <- import(bg, format = "bedGraph", seqinfo = seqinfo(txdb))
  strand(gr) <- ifelse(score(gr) >= 0, "+", "-")
  score(gr) <- abs(score(gr))
  gr_fw <- sort(gr[strand(gr) == "+"])
  gr_rev <- sort(gr[strand(gr) == "-"])
  export(gr_fw, sub(".bedgraph.gz", "_Plus.bw", bg), format = "BigWig")
  export(gr_rev, sub(".bedgraph.gz", "_Minus.bw", bg), format = "BigWig")
}

# Draw scatterplot of reproducibility between biological replicates (Fig. S1D):
all_tiles <- tileGenome(seqinfo(txdb), tilewidth = 10, cut.last.tile.in.chrom = TRUE)

rep1 <- list.files(bamdir, pattern = "biorep1.*bedgraph.gz$")
rep2 <- sub("biorep1", "biorep2", rep1)
data_rep1 <- batchReadTrackData(rep1, dir = bamdir, format = "bedGraph", seqinfo  = seqinfo(txdb))
data_rep2 <- batchReadTrackData(rep2, dir = bamdir, format = "bedGraph", seqinfo  = seqinfo(txdb))

genes <- genes(txdb, columns = c("gene_id", "tx_type"))
mcols(genes)$tx_type <- unname(unlist(mcols(genes)$tx_type))
genes_npcd <- genes[mcols(genes)$tx_type == "protein_coding" & seqnames(genes) %in% 1:5] # nuclear protein-coding genes in Plantsmart28
genes_npcd_ext <- suppressWarnings(trim(resize(genes_npcd, width(genes_npcd) + 100, "end")))
genes_npcd_ext <- suppressWarnings(trim(resize(genes_npcd_ext, width(genes_npcd_ext) + 500, "start")))
genes_npcd_m <- GenomicRanges::reduce(genes_npcd_ext)

data_rep1 <- endoapply(data_rep1, normalizeGR, by = genes_npcd_m) # normalize track to 1M tags in extended nuclear protein-coding genes
data_rep2 <- endoapply(data_rep2, normalizeGR, by = genes_npcd_m)

for (i in seq_along(rep1)) {
  data <- list(data_rep1[[i]], data_rep2[[i]])
  names(data) <- c(rep1[[i]], rep2[[i]])
  mat <- getOverlappingScores(all_tiles, data, value = "count_matrix")
  mat <- mat[mat[, 1] > 0 & mat[, 2] > 0, ]
  df <- as.data.frame(mat)
  s1 <- names(df)[[1]]
  s2 <- names(df)[[2]]
  r <- round(cor(df[, 1], df[, 2], method = "pearson"), 3)
  ttl <- paste0("Technical reproducibility ", sub(".bedgraph.gz", "", s1), " vs ", sub(".bedgraph.gz", "", s2), " (10 bp windows, all genes)")
  p <- ggplot(df, aes(x = eval(parse(text = s1)), y = eval(parse(text = s2)))) + geom_point(shape = 19, size = 1.5, alpha = 0.01) + 
    ggtitle(paste0("Pearson r = ", r)) + xlab(s1) + ylab(s2) + xlim(0, 50) + ylim(0, 50)
  for (ext in c(".pdf", ".png")) { ggsave(paste0(ttl, ext), plot = p, width = 7, height = 7, units = "in") }
}