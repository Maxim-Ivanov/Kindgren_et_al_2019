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

scripts <- c("convertGAlignmentsToCoverage.R", "saveGRangesAsBedGraph.R", "mergeBedgraphs.R", "loadNETSeqBAM.R", "removeFirstAndLastExons.R")
for (script in scripts) { source(script) }

findSpliceSites <- function(ebg) {
  res <- removeFirstAndLastExons(ebg)
  exons_wo_first <- res[[1]]
  exons_wo_last <- res[[2]]
  ss_5p <- reduce(resize(resize(exons_wo_last, 1, "end"), 3, "center")) # the last base of exon (+ 1 base up  + 1 base down)
  ss_3p <- reduce(resize(flank(exons_wo_first, 1), 3, "center")) # the last base of the upstream intron (+ 1 base up  + 1 base down)
  return(list(ss_5p, ss_3p))
}

ss1 <- findSpliceSites(ebg1)
donor_1 <- ss1[[1]]
acc_1 <- ss1[[2]]
ss2 <- findSpliceSites(ebg3)
donor_2 <- ss2[[1]]
acc_2 <- ss2[[2]]
all_donor <- sort(reduce(c(donor_1, donor_2))) # all donor splice sites
all_acc <- sort(reduce(c(acc_1, acc_2))) # all acceptor splice sites

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
