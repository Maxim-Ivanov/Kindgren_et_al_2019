# The script below was used to compute both the promoter-proximal stalling index (TSS SI) and the intronic stalling index (ISI);
# This file also contains the code to reproduce figures 8b and S6c;

##### PART 1: promoter-proximal stalling index (TSS SI) #####

# The methodology for calculating promoter-proximal RNAPII stalling index was inspired by Zhu 2018 (PMID 30374093):
# We expect to observe the TSS stalling event somewhere in the interval [TSS-100, TSS+300];
# A sliding window of fixed width (100 bp) moves in 10 bp steps along this interval;
# PlaNET-Seq coverage (averaged between all samples) is calculated for each position of the window;
# Position with the strongest signal is considered as the most representative position for TSS stalling in given gene;
# (if multiple positions had the same highest signal, then the position closest to [TSS+100] was chosen);
# Finally, stalling index is calculated as the plaNET-Seq FPKM at the TSS stalling position divided by FPKM in the coding region [TSS+300, PAS-100];

library(SummarizedExperiment)
library(GenomicFeatures)
library(edgeR)
library(rtracklayer)
library(ggplot2)
library(dplyr)
library(reshape2)

r_dir <- "." # change to the folder where you saved the custom functions from https://github.com/Maxim-Ivanov/Kindgren_et_al_2019
scripts <- c("getOverlappingScores.R", "windowsFromGRanges.R", "findBestWindow.R", "normalizeGR.R")
for (script in scripts) { source(file.path(r_dir, script)) }

genes_araport_adj <- readRDS("genes_araport_adj.RDS") # see 04-Adjustment_Araport11.R
genes_npcd <- genes_araport_adj[mcols(genes_araport_adj)$tx_type == "mRNA" & seqnames(genes_araport_adj) %in% 1:5]

txdb_araport <- makeTxDbFromGFF("Araport11_GFF3_genes_transposons.201606.gff.gz") # https://www.arabidopsis.org/download_files/Genes/Araport11_genome_release/Araport11_GFF3_genes_transposons.201606.gff.gz
ebg_araport <- exonsBy(txdb_araport, by = "gene") # exons grouped by gene
seqinfo(ebg_araport, new2old = as.integer(c(1:5, 7, 6))) <- seqinfo(genes_araport_adj)

library(BSgenome.Athaliana.TAIR.TAIR9)
bsgen <- BSgenome.Athaliana.TAIR.TAIR9
seqlevels(bsgen) <- seqlevels(genes_araport_adj)

# Load plaNET-Seq DE data:
genes_de <- readRDS("genes_de.RDS") # see 07-DESeq2_pipeline.R
genes_de_gr <- rowRanges(genes3_de)
genes_cm <- assay(genes3_de)

# Load plaNET-Seq data (normalized to 1M tags in nuclear protein-coding genes):
planet_data <- readRDS("PlaNET-Seq_data_norm1M.RDS") # see 05-Metagenes.R

# Calculate FPKM transcription of whole genes in the wild type sample:
cov <- as.numeric(getOverlappingScores(genes_de_gr, planet_data[1], value = "count_matrix"))
mcols(genes_de_gr)$fpkm_wt <- round(cov / width(genes_de_gr) * 1000, 5)

# Subset Araport11 genes to width >= 1 Kb, WT FPKM >= 1 and no overlap with other genes:
genes_ext <- suppressWarnings(trim(resize(granges(genes_de_gr), width(genes_de_gr) + 100, "end")))
good <- width(genes_de_gr) >= 1000 & mcols(genes_de_gr)$fpkm_wt >= 1 & countOverlaps(genes_ext, genes_de_gr) == 1
genes_good <- genes_de_gr[good] # n = 17195

# Find the position of TSS stalling peak:
prom <- trim(resize(resize(granges(genes_good), 300, "start"), 400, "end"))
win <- windowsFromGRanges(prom, window_width = 100, window_offset = 10)
wcov <- getOverlappingScores(win, planet_data, value = "count_matrix")
mcols(win)$cov <- round(rowMeans(wcov), 5) # average plaNET-Seq signal among all samples
win_c <- resize(win, 1, "center")
prior <- resize(resize(granges(genes_good), 100, "start"), 1, "end") # [TSS+100] was taken as the best prior estimate for the TSS stalling position
prior_par <- prior[mcols(win)$ID]
dist <- start(win_c) - start(prior_par)
mcols(win)$dist <- ifelse(strand(win_c) == "+", dist, -dist) # negative distance: window is upstream from the gene TSS
mcols(win)$idx <- 1:length(win)
best_win_idx <- by(mcols(win)[, c("cov", "dist", "idx")], INDICES = mcols(win)$ID, FUN = findBestWindow, simplify = FALSE) %>% unlist()
w1 <- win[best_win_idx]

# PlaNET-Seq coverage of the TSS stalling peak:
cov1 <- getOverlappingScores(w1, planet_data, value = "count_matrix")

# PlaNET-Seq coverage in coding regions of respective genes:
w2 <- resize(granges(genes_good), width(genes_good) - 300, "end")
w2 <- resize(w2, width(w2) - 100, "start")
cov2 <- getOverlappingScores(w2, planet_data, value = "count_matrix")
cov2_norm <- cov2 / width(w2) * 100 # normalize coverage to 100 bp of gene width

# Compute TSS stalling index:
tss_si <- cov1 / cov2_norm
tss_si[is.nan(tss_si)] <- 0
tss_si[is.infinite(tss_si)] <- 0
colnames(tss_si) <- sub("plaNET", "TSS_SI", colnames(tss_si))
tss_si <- round(tss_si, 3)

# Combine TSS SI results with the DE data:
mc <- mcols(genes_good)
mc_de <- mc[, grepl("_de", colnames(mc))]
df <- as.data.frame(cbind(mc_de, tss_pi))

for (i in grep("_de", colnames(df))) {
  df[, i] <- factor(df[, i], levels = c("No", "Up", "Down"))
}

# Fig. 8b (boxplot of WT TSS SI in genes which are Up/Down/nonDE in Cold 3h):
pi_col <- "TSS_SI_WT"
de_col <- "Cold_3h_vs_WT_de"
ttl <- "TSS Stalling Index in WT vs DE status in Cold 3h (Fig. 8b)"
stats <- tapply(df[, pi_col], df[, de_col], boxplot.stats)
max_y <- max(unlist(lapply(stats, function(x) { return(x[[1]][[5]]) })))
med <- median(df[df[, de_col] == "No", pi_col])
p <- ggplot(df, aes(x = eval(parse(text = de_col)), y = eval(parse(text = pi_col)))) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0, max_y)) + 
  ggtitle(ttl) + geom_hline(yintercept = med, linetype = "dashed")
for (ext in c(".png", ".pdf")) { ggsave(paste0(ttl, ext), plot = p, width = 7, height = 7, units = "in") }

# Number of genes in each group shown on the boxplot:
table(df[, de_col])

# Statistical significance of the TSS SI differences between groups:
x <- df[df[, de_col] == "No", pi_col]
y1 <- df[df[, de_col] == "Up", pi_col]
y2 <- df[df[, de_col] == "Down", pi_col]
wilcox.test(x, y1)$p.value # upregulated genes vs non-regulated ones
wilcox.test(x, y2)$p.value # downregulated genes vs non-regulated ones


##### PART 2: intronic stalling index (ISI) #####
# Step 1: find positions of intronic RNAPII stalling in each intron (in nuclear protein-coding genes) using pNET-Seq Ser5P sample;
# This sample was chosen because it was demonstrated to have the most profound intronic peak among all plaNET-Seq and pNET-Seq samples (see Fig. 4f-g, lower panels);
# Position of the "best" window is found using the same approach as above for TSS SI;
# Step 2: calculate intron stalling index in different plaNET samples using the coordinates of "best" windows found on step 1;

ww <- 10 # width of the sliding window (bp)
step <- 5 # step for sliding (5 bp)

# Load pNET-Seq data:
pnet_dir <- "." # change to the directory containing merged pNET-Seq Bedgraph files (obtained from 02-Postprocessing_plaNET-Seq.R)
pnet_files <- list.files(pnet_dir, pattern = "merged_fw_rev.bedgraph.gz$")
pnet_data <- batchReadTrackData(pnet_files, dir = pnet_dir, format = "bedGraph", seqinfo = seqinfo(genes_araport_adj))
names(pnet_data) <- paste0("pNET_", sub("_merged_fw_rev.bedgraph.gz", "", names(pnet_data)))
ser5p <- pnet_data[[2]] # raw pNET-Seq tags
ser5p_norm <- normalizeGR(ser5p) # pNET-Seq signal normalized to 1M tags

# Extract and enumerate introns:
ibg <- psetdiff(unlist(range(ebg)), ebg) # introns grouped by gene
introns <- unlist(ibg) # n = 119681
inum <- lengths(ibg) # number of introns in each gene
inum_nonzero <- inum[inum > 0] # exclude intronless genes
ibg_enum_fw <- as.vector(unlist(lapply(inum_nonzero, function(x) { seq(1, x) })))
ibg_enum_rev <- as.vector(unlist(lapply(inum_nonzero, function(x) { seq(x, 1) })))
mcols(introns)$pos_fw <- ifelse(strand(introns) == "+", ibg_enum_fw, ibg_enum_rev) # intron number from gene start
mcols(introns)$pos_rev <- ifelse(strand(introns) == "+", -ibg_enum_rev, -ibg_enum_fw) # intron number from gene end
mcols(introns)$inum <- rep(inum_nonzero, inum_nonzero) # number of introns in given gene
mcols(introns)$gene <- rep(unlist(range(ebg))[inum > 0], inum_nonzero) # gene coordinates

# Calculate FPKM values for both introns and thei host genes in the pNET-Seq Ser5P sample:
fpkm_genes <- round(as.numeric(getOverlappingScores(genes_araport_adj, list("Ser5P" = ser5p_norm), value = "count_matrix")) / width(genes_araport_adj) * 1000, 3)
mcols(introns)$gene_fpkm <- rep(fpkm_genes[inum > 0], inum_nonzero)
mcols(introns)$intron_fpkm <- round(as.numeric(getOverlappingScores(introns, list("Ser5P" = ser5p_norm), value = "count_matrix")) / width(introns) * 1000, 3)

# Choose introns to calculate the ISI:
introns_good <- introns[width(introns) >= 50 & width(introns) <= 300 & countOverlaps(introns, genes_araport_adj) == 1 & names(introns) %in% names(genes_npcd)] # n = 102774

# Find the total number of raw tags on each intron:
mcols(introns_good)$intron_tags <- as.integer(getOverlappingScores(introns_good, list("Ser5P" = ser5p), value = "count_matrix"))

# Find the "best" window for each intron (using pNET-Seq Ser5P data):
win <- windowsFromGRanges(introns_good, window_width = ww, window_offset = step)
mcols(win)$cov <- as.integer(getOverlappingScores(win, list("Ser5P" = ser5p), value = "count_matrix"))
mcols(win)$idx <- 1:length(win)
best_win_idx <- unlist(by(mcols(win)[, c("cov", "idx")], INDICES = mcols(win)$ID, FUN = findBestWindow, random = TRUE, simplify = FALSE))
win_best <- win[best_win_idx]

# Calculate ISI:
mcols(introns_good)$win_best <- granges(win_best)
mcols(introns_good)$win_best_tags <- mcols(win_best)$cov
isi <- round(mcols(win_best)$cov / mcols(introns_good)$intron_tags * (width(introns_good) / width(win_best)), 3)
isi[is.nan(isi)] <- 0
mcols(introns_good)$ISI <- isi

# Remove low transcribed introns (which are probable to have an artifactually overestimated ISI):
introns_good <- introns_good[mcols(introns_good)$intron_fpkm >= 10] # n = 14795

# Stratify the remaining introns by ISI into "strong", "medium" and "weak":
strong <- mcols(introns_good)$ISI >= 5.5 # n = 3704
weak <- mcols(introns_good)$ISI <= 3.5 # n = 4328
medium <- !strong & !weak # n = 6763
groups <- rep("Medium", length(introns_good))
groups[strong] <- "Strong"
groups[weak] <- "Weak"
groups <- factor(groups, levels = c("Weak", "Medium", "Strong"))
mcols(introns_good)$groups <- groups

# Fig. S6c (boxplot of intron lengths grouped by the ISI):
df <- data.frame("Width" = width(introns_good), "Group" = groups)
ttl <- "Intron width vs ISI (Fig. S6c)"
p <- ggplot(df, aes(x = Group, y = Width, fill = Group)) + geom_boxplot() + ggtitle(ttl)
for (ext in c(".pdf", ".png")) { suppressWarnings(ggsave(paste0(ttl, ext), plot = p, width = 7, height = 7, units = "in")) }
