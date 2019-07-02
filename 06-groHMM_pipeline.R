# This script describes how transcripts were called de novo from the plaNET-Seq data;
# It also allows to reproduce Figures 5B, 5D, 5E, 6B, 7B, 7D, S6B, S6C, S6G; 

library(groHMM)
options(mc.cores=getCores(4))
library(dplyr)
library(rtracklayer)
library(ggplot2)

r_dir <- "." # change to the folder where you saved the custom functions from https://github.com/Maxim-Ivanov/Kindgren_et_al_2019
scripts <- c("batchReadTrackData.R", "getOverlappingScores.R", "poverlapsGRanges.R", "parallelOverlapType.R", "normalizeGR.R", "annotateNoncodingTranscripts.R")
for (script in scripts) { source(file.path(r_dir, script)) }

# Load adjusted Araport11 genes (see 04-Adjustment_Araport11.R for details):
genes_araport_adj <- readRDS("genes_araport_adj.RDS")

# Flip genes to the opposite strand:
genes_AS <- genes_araport_adj
strand(genes_AS) <- ifelse(strand(genes_AS) == "+", "-", "+")

# Load the original Araport11 annotation:
txdb_araport <- makeTxDbFromGFF("Araport11_GFF3_genes_transposons.201606.gff.gz")
ebg_araport <- exonsBy(txdb_araport, by = "gene") # exons grouped by gene
seqinfo(ebg_araport, new2old = as.integer(c(1:5, 7, 6))) <- seqinfo(genes_araport_adj)

# Extract protein-coding genes (for normalization of plaNET-Seq tracks):
genes_npcd <- genes_araport_adj[mcols(genes_araport_adj)$tx_type == "mRNA" & seqnames(genes_araport_adj) %in% 1:5] # nuclear protein-coding genes
genes_npcd_ext <- suppressWarnings(trim(resize(genes_npcd, width(genes_npcd) + 100, "end"))) # extend by 100 bp upstream
genes_npcd_ext <- suppressWarnings(trim(resize(genes_npcd_ext, width(genes_npcd_ext) + 500, "start")))  # extend by 500 bp downstream to capture pA peaks in plaNET-Seq data
genes_npcd_m <- reduce(genes_npcd_ext) # merge overlapping intervals

# Load plaNET-Seq files (merged biological replicates, raw tag counts):
planet_dir <- "." # change to the directory containing merged plaNET-Seq Bedgraph files obtained from 02-Postprocessing_plaNET-Seq.R
planet_files <- list.files(planet_dir, pattern = "merged_fw_rev.bedgraph.gz$")
planet_data <- batchReadTrackData(planet_files, dir = planet_dir, format = "bedGraph", seqinfo = seqinfo(genes_araport_adj))
names(planet_data) <- paste0("plaNET_", sub("_merged_fw_rev.bedgraph.gz", "", names(planet_data)))
planet_data <- lapply(planet_data, function(gr) { seqlevels(gr) <- seqlevelsInUse(gr); return(gr) }) # due to unknown reason, groHMM does not like unused seqlevels...

# Normalize to 1M reads in nuclear protein-coding genes (for coverage calculation only, not for groHMM):
planet_data_norm1M <- endoapply(planet_data, normalizeGR, by = genes_npcd_m)

# Run GroHMM on every plaNET-Seq genotype/condition:
planet_grohmm <- vector("list", length(planet_data))
names(planet_grohmm) <- names(planet_data)
for (i in seq_along(planet_data)) {
  data <- planet_data[[i]]
  # Custom groHMM settings for the gene-dense Arabidopsis genome:
  tx <- detectTranscripts(data, LtProbB = -15, UTS = 5, threshold = 1, mc.silent = TRUE)$transcripts
  mcols(tx) <- NULL
  mcols(tx)$sample <- i
  tx <- tx[width(tx) >= 200]
  planet_grohmm[[i]] <- tx
}

# Export transcripts called in each sample as BED files:
for (i in seq_along(planet_grohmm)) {
  export(planet_grohmm[[i]], paste0("GroHMM_transcripts_", names(planet_grohmm)[[i]], ".bed"), format = "BED")
}

# Cluster and merge strongly overlapping intervals which were called in different conditions:
findUniqueBestOverlap <- function(x) {
  one_best <- vector("logical", length(x))
  one_best[which.max(x)] <- TRUE
  return(one_best)
}

clusterByReciprocalOverlap <- function(x, y, threshold = 0.5) {
  x <- sort(x)
  y <- sort(y)
  if (!"cluster" %in% names(mcols(x))) {
    mcols(x)$cluster <- seq(1, length(x))
  }
  hits <- findOverlaps(y, x)
  x2 <- x[subjectHits(hits)]
  y2 <- y[queryHits(hits)]
  over <- width(pintersect(x2, y2)) / width(punion(x2, y2))
  best_over <- unlist(tapply(over, list(queryHits(hits)), findUniqueBestOverlap))
  valid_over <- over >= threshold
  best_valid <- best_over & valid_over
  x3 <- x2[best_valid]
  y3 <- y2[best_valid]
  y_only <- y[!(y %in% y3)]
  mcols(y3)$cluster <- mcols(x3)$cluster
  mcols(y_only)$cluster <- seq(max(mcols(x)$cluster) + 1, length.out = length(y_only))
  out <- c(x, y3, y_only)
  return(out)
}

mergeClusters <- function(gr) {
  gr <- gr[order(mcols(gr)$cluster)] # sort by cluster number
  grl <- split(gr, mcols(gr)$cluster)
  n_samples = max(mcols(gr)$sample)
  mat <- matrix(0, nrow = length(grl), ncol = n_samples)
  chroms <- unlist(runValue(seqnames(grl)))
  strands <- unlist(runValue(strand(grl)))
  starts <- unlist(lapply(start(grl), median))
  ends <- unlist(lapply(end(grl), median))
  smp <- tapply(mcols(gr)$sample, list(mcols(gr)$cluster), c)
  for (i in 1:nrow(mat)) {
    mat[i, smp[[i]]] <- 1
  }
  df <- as.data.frame(mat)
  names(df) <- paste0("s", 1:n_samples)
  res <- GRanges(seqnames = chroms, ranges = IRanges(starts, end = ends), strand = strands)
  mcols(res) <- df # shows in which plaNET-Seq samples the transcript was called (encoded by 1 or 0)
  return(sort(res))
}

clustered <- Reduce(clusterByReciprocalOverlap, planet_grohmm)
merged <- mergeClusters(clustered)
seqinfo(merged) <- seqinfo(genes_araport_adj)

# Export merged groHMM transcripts as BED file:
export(merged, "GroHMM_transcripts_merged_across_samples.bed", format = "BED")
saveRDS(merged, "GroHMM_transcripts_merged_across_samples.RDS")

# Calculate some useful statistics:
# 1) Number of samples where the transcript was called ("positive" samples):
posSamples <- rowSums(as.matrix(mcols(merged)))
# 2) plaNET-Seq FPKM of transcripts in each sample:
fpkm <- round(getOverlappingScores(merged, planet_data_norm1M, value = "count_matrix") / width(merged) * 1000, 5)
colnames(fpkm) <- paste0("FPKM_", colnames(fpkm))
# 3) Mean FPKM in "positive" samples:
score <- rowMeans(fpkm[as.logical(mcols(merged))])
# Add the results to mcols(merged):
mcols(merged) <- cbind(mcols(merged), data.frame("posSamples" = possamples), as.data.frame(fpkm), data.frame("score" = score))


##### ANNOTATE GROHMM TRANSCRIPTS BY OVERLAP WITH KNOWN GENES IN ARAPORT11 #####

# Count overlaps with genes on the same strand:
over_S <- countOverlaps(merged, genes_araport_adj)
tx_over_S <- merged[over_S > 0]
tx_noOver_S <- merged[over_S == 0]
# Prepare the output data frame:
df <- data.frame("over_S" = rep(FALSE, length(merged)))
df$over_S[over_S > 0] <- TRUE # TRUE if transcript overlaps any gene on the same strand
df$N_over_S <- over_S # number of overlaps with genes on the same strand

# For transcripts overlapping at least one gene:
# 1) Decide if the overlap is strong (reciprocal 50%):
hits <- findOverlaps(tx_over_S, genes_araport_adj)
tx_over_S_par <- tx_over_S[queryHits(hits)]
genes_par <- genes_araport_adj[subjectHits(hits)]
overlap <- width(pintersect(tx_over_S_par, genes_par))
overlap_tx <- overlap / width(tx_over_S_par)
overlap_genes <- overlap / width(genes_par)
valid_overlap <- overlap_tx >= 0.5 | overlap_genes >= 0.5 # either transcript or gene is covered by at least 50% of its width
strong_over <- as.logical(unlist(tapply(valid_overlap, queryHits(hits), any)))
df$strong_over_S <- FALSE
df$strong_over_S[df$over_S] <- strong_over # TRUE if transcript has at least one strong overlap on the same strand
# 2) Detect upstream and downstream overlaps:
a <- start(tx_over_S_par) <= start(genes_par) & end(tx_over_S_par) < end(genes_par)
c <- start(tx_over_S_par) > start(genes_par) & end(tx_over_S_par) >= end(genes_par)
up <- ifelse(strand(tx_over_S_par) == "+", a, c)
over_up <- as.logical(unlist(tapply(up, queryHits(hits), any)))
df$over_S_up <- FALSE
df$over_S_up[df$over_S] <- over_up # TRUE if transcript has a 5' overlap on the same strand
down <- ifelse(strand(tx_over_S_par) == "+", c, a)
over_down <- as.logical(unlist(tapply(down, queryHits(hits), any)))
df$over_S_down <- FALSE
df$over_S_down[df$over_S] <- over_down # TRUE if transcript has a 3' overlap on the same strand

# For transcripts not overlapping genes on the same strand:
# 1) Get the distance to the nearest upstream gene:
idx_up <- follow(tx_noOver_S, genes_araport_adj) # can contain NA
nearest_up <- genes_araport_adj[idx_up[!is.na(idx_up)]]
par_up <- tx_noOver_S[!is.na(idx_up)]
dist_up <- width(pgap(nearest_up, par_up))
dist_up_full <- rep(NA, length(tx_noOver_S)) # expand to the original length (with NA)
dist_up_full[!is.na(idx_up)] <- dist_up
df$dist_up <- NA
df$dist_up[!df$over_S] <- dist_up_full # distance to the nearest upstream gene on the same strand
# 2) Do the same for the nearest downstream gene:
idx_down <- precede(tx_noOver_S, genes_araport_adj)
nearest_down <- genes_araport_adj[idx_down[!is.na(idx_down)]]
par_down <- tx_noOver_S[!is.na(idx_down)]
dist_down <- width(pgap(nearest_down, par_down))
dist_down_full <- rep(NA, length(tx_noOver_S))
dist_down_full[!is.na(idx_down)] <- dist_down
df$dist_down <- NA
df$dist_down[!df$over_S] <- dist_down_full # distance to the nearest downstream gene on the same strand

# Count overlaps with genes on the opposite strand:
over_AS <- countOverlaps(merged, genes_AS)
df$over_AS <- FALSE
df$over_AS[over_AS > 0] <- TRUE # TRUE if transcript overlaps any gene on the opposite strand
df$N_over_AS <- over_AS # number of overlaps with genes on the opposite strand

# Detect DNC transcripts (start within 500 bp from TSS of a known gene):
dnc_intervals <- suppressWarnings(trim(flank(genes_AS, 500, start = FALSE)))
tx_starts <- resize(merged, 1, fix = "start")
dnc <- tx_starts %over% dnc_intervals
df$DNC <- dnc # TRUE is transcript is divergent with respect to a known gene on the opposite strand

# Detect transcripts which start inside of a gene in the antisense orientation:
inside <- countOverlaps(tx_starts, genes_AS)
df$start_within_AS <- inside # TRUE if transcript starts inside of a known gene on the opposite strand
# For these transcripts, calculate the distance from gene TSS to the transcript start:
tx_inside <- merged[inside > 0]
hits2 <- findOverlaps(tx_inside, genes_AS)
tx_inside_par <- tx_inside[queryHits(hits2)]
genes_AS_par <- genes_AS[subjectHits(hits2)]
tss_AS_par <- resize(genes_AS_par, 1, "end")
tx_inside_start_par <- resize(tx_inside_par, 1, "start")
gaps_AS <- pgap(tx_inside_start_par, tss_AS_par)
dist_AS <- width(gaps_AS)
# If there are multiple known genes per given DNC transcript, choose the one with the minimal gap:
min_dist <- unlist(tapply(dist_AS, queryHits(hits2), function(x) { out <- vector("logical", length(x)); out[which.min(x)] <- TRUE; return(out) }))
min_gaps_AS <- gaps_AS[min_dist]
min_dist_AS <- dist_AS[min_dist]
df$dist_TSS_start_AS <- NA
df$dist_TSS_start_AS[df$start_within_AS > 0] <- min_dist_AS # distance between divTSS and coding TSS
# Check if the gap between TSS of the known gene and the start of transcript on the antisense strand overlaps with any nucleosome:
nps <- read.table("Ath_leaf_NPS.gff", sep = "\t", header = FALSE) # download from http://plantdhs.org/static/download/Ath_leaf_NPS.gff.gz
nps <- GRanges(seqnames=sub("Chr", "", nps$V1), ranges=IRanges(nps$V4, end=nps$V5), seqinfo=seqinfo(genes_araport_adj))
nps_center <- resize(nps, 1, "center")
over_nucl <- countOverlaps(min_gaps_AS, nps_center)
df$over_nucl_AS <- NA
df$over_nucl_AS[df$start_within_AS > 0] <- over_nucl # number of nucleosomes between divTSS and coding TSS

# Classify transcripts based on the collected values:
#df$conv <- df$start_within_AS > 0 & df$dist_TSS_start_AS <= 500 & df$over_nucl_AS == 0
#df$intragenic_AS <- df$start_within_AS > 0 & !df$conv
#df$outer_AS <- isTRUE(over_AS) & !df$conv & !df$intragenic_AS

# Add the data frame to mcols(merged):
mcols(merged) <- cbind(mcols(merged), df)

# Remove transcripts with strong overlap to known genes:
novel <- merged[!mcols(merged)$strong_over_S] # n = 8954
# (the transcrripts retained are considered novel)

# Mark transcripts which are subintervals of longer intervals:
skipSubintervals <- function(gr, threshold = 0.9) {
  hits <- findOverlaps(gr, gr)
  hits <- hits[queryHits(hits) != subjectHits(hits)]
  x <- gr[queryHits(hits)]
  y <- gr[subjectHits(hits)]
  xy <- width(pintersect(x, y))
  bad <- xy / width(x) >= threshold
  # Interval is skipped if it has at least one really strong overlap (>= 90%) with another interval:
  skip <- as.logical(tapply(bad, list(queryHits(hits)), any))
  gr_noOver <- gr[!gr %in% x]
  gr_over <- gr[gr %in% x]
  gr_good <- gr_over[!skip]
  subintervals <- !(gr %in% c(gr_noOver, gr_good)) # TRUE for subintervals
  return(subintervals)
}

mcols(novel)$Subint <- skipSubintervals(novel) # FALSE for intervals which are not subintervals (n = 7228)

# Export novel transcripts as BED file:
export(novel, "Novel_groHMM_transcripts.bed", format = "BED")
saveRDS(novel, "Novel_groHMM_transcripts.RDS")

##### Fig. 5B: classification of novel transcripts and known genes #####

# Classify each transcript as "Genic" (strong overlap with any known gene on the same strand), 
# "Upstream" (weak overlap with 5' end of a known gene),
# "Downstream" (weak overlap with 3' end),
# "Divergent" (transcript starts within 500 bp from known gene on the opposite strand),
# "Convergent" (transcript overlaps a known gene on the opposite strand and starts within the first 50% of the gene length),
# "TTS_AS" (the same as convergent but starts within 50%-120% of gene length),
# "Distal_AS" (the same as TSS_AS but starts further downstream of the known gene),
# "Intergenic" (transcripts without overlap with any known gene on both strands);

# Classify groHMM transcripts:
mat1 <- annotate_noncoding_transcripts(merged, genes_araport_adj) # the logical matrix contains TRUE/FALSE for each transcript in each category
# Skip genic transcripts:
non_genic <- !mat1[, "Genic"]
mat1 <- mat1[non_genic, ] # n = 8954
# Skip subintervals:
mat1 <- mat1[!skipSubintervals(merged[non_genic]), ] # n = 7228

# Classify known non-coding transcripts in Araport11:
pc <- mcols(genes_araport_adj)$tx_type == "mRNA" # protein-coding
snc <- mcols(genes_araport_adj)$tx_type %in% c("tRNA", "rRNA", "snRNA", "snoRNA") # short non-coding
lnc <- !pc & !snc # long non-coding
genes_lnc <- genes_araport_adj[lnc]
genes_wo_lnc <- genes_araport_adj[!lnc]
mat2 <- annotate_noncoding_transcripts(genes_lnc, genes_wo_lnc) # n = 5565
# For consistency, skip genic intervals:
mat2 <- mat2[!mat2[, "Genic"], ] # n = 5412

# The matrix can contain multiple annotations for the same gene. The only non-conflicting annotation is Intergenic (it is TRUE only when all others are FALSE).
# One of the possible strategies is to invent an arbitrary hierarchy of annotations.
# If a gene belongs to multiple categories at once, the leftmost category is used:
hierarchy <- c("Genic", "Divergent", "Convergent", "TTS_AS", "Distal_AS", "Up", "Down", "Intergenic") 

apply_hierarchy <- function(mat, hier) {
  out <- rep(NA, nrow(mat))
  for (i in seq_along(hier)) {
    name <- hier[[i]]
    out[mat[, name] & is.na(out)] <- name
  }
  return(out)
}


h1 <- apply_hierarchy(mat1, hierarchy)
df1 <- as.data.frame(table(h1))
names(df1) <- c("Type", "Count")
df1$Group <- "Novel"
h2 <- apply_hierarchy(mat2, hierarchy)
df2 <- as.data.frame(table(h2))
names(df2) <- c("Type", "Count")
df2$Group <- "Known"

df3 <- rbind(df1, df2)
df3 <- subset(df3, Type != "Genic")
df3$Type <- factor(df3$Type, levels = c("Up", "Down", "Divergent", "Convergent", "TTS_AS", "Distal_AS", "Intergenic"))
df3$Group <- factor(df3$Group, levels = c("Novel", "Known"))
p2 <- ggplot(df3, aes(x = Type, y = Count, fill = Group)) + geom_bar(stat = "identity") + geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 3)
for (ext in c(".pdf", ".png")) { ggsave(paste0("Novel and known ncRNA classification", ext), plot = 2, width = 7, height = 7, units = "in") }


##### ANALYSIS OF DIVERGENT NOVEL TRANSCRIPTS #####

# Extract divergent transcripts:
dnc <- subset(novel, DNC & !Subint) # n = 917
# Find their mate genes:
idx <- follow(dnc, genes_AS)
mate_genes <- genes_araport_adj[idx] 
# Subset to nuclear protein-coding mate genes:
npcd <- mcols(mate_genes)$tx_type == "mRNA" & seqnames(mate_genes) %in% 1:5
dnc <- dnc[npcd]
mate_genes <- mate_genes[npcd]
# Deduplicate on mate genes:
dupl <- duplicated(names(mate_genes))
dnc <- dnc[!dupl]
mate_genes <- mate_genes[!dupl]

# Fig. 5D: histogram/density of abs. distances between divTSS and coding TSS:
dist <- width(pgap(dnc, mate_genes3, ignore.strand = TRUE))
df1 <- data.frame("Dist" = dist)
p1 <- ggplot(df1, aes(x = Dist, y = ..density..)) + geom_histogram(bins = 30) + geom_density(colour = "red") + xlab("Distance, bp") + ggtitle("Absolute distance to TSS of the mate gene")
for (ext in c(".png", ".pdf")) { ggsave(paste0("Distance to the mate gene", ext), plot = p1, width = 6, height = 6) }

# Fig. 6B: Boxplot of plaNET-Seq FPKM (WT sample) in DNC vs non-DNC genes:
genes_npcd_fpkm <- as.numeric(getOverlappingScores(genes_npcd, planet_data[1], value = "count_matrix")) / width(genes_npcd) * 1000
silent <- genes_npcd_fpkm < 1
df2 <- data.frame("Expr" = genes_npcd_fpkm, "Group" = ifelse(names(genes_npcd) %in% names(mate_genes), "With_DNC", "Without_DNC"), "Silent" = silent)
df2_expr <- subset(df2, !Silent)
df2_expr$Expr <- log2(df2_expr$Expr)
p2 <- ggplot(df2_expr, aes(x = Group, y = Expr)) + geom_boxplot() + ylab("Log2 FPKM") + ggtitle("Expressed genes (plaNET-Seq WT FPKM >= 1)")
for (ext in c(".png", ".pdf")) { ggsave(paste0("Expression of genes with vs without DNC", ext), plot = p2, width = 6, height = 6) }
# Mann-Whitney p-value:
wilcox.test(x = df2_expr[df2_expr$Group == "With_DNC", "Expr"], y = df2_expr[df2_expr$Group == "Without_DNC", "Expr"])$p.value # 5.798428e-23
sum(df2_expr$Group == "With_DNC") # n = 711
sum(df2_expr$Group == "Without_DNC") # n = 19300

##### ANALYSIS OF CONVERGENT NOVEL TRANSCRIPTS #####

# Extract all AS trancripts (either start within the host gene, or at least overlap its 3' end):
as <- subset(novel, over_AS & N_over_AS == 1 & !Subint) # n = 5313
# Flip them to the opposite (coding) strand:
as_flipped <- as
strand(as_flipped) <- ifelse(strand(as) == "+", "-", "+")
as_flipped_tss <- resize(as_flipped, width = 1, fix = "end")
# Find their host genes:
hits <- findOverlaps(as_flipped, genes_araport_adj, select = "first")
host_genes <- genes_araport_adj[hits]
# Absolute and relative positions of AS transcripts within their host genes:
abs_pos <- width(pgap(resize(host_genes, 1, "start"), as_flipped_tss))
rel_pos <- abs_pos / width(host_genes)
mcols(as)$rel_pos <- rel_pos # can exceed 1.0
# Skip transcripts with relative distance > 1.2:
as <- as[mcols(as)$rel_pos <= 1.2]
host_genes <- host_genes[mcols(as)$rel_pos <= 1.2]
# Classify AS transcripts as "convergent" and "TSS-AS":
mcols(as)$category <- ifelse(mcols(as)$rel_pos <= 0.5, "Convergent", "TSS_AS")

### Fig. 5E (histogram of relative distances between AS TSS and coding TSS):
df1 <- data.frame("Relative_position" = mcols(as)$rel_pos)
p <- ggplot(df1, aes(x = Relative_position)) + geom_histogram(breaks = seq(0, 1.2, by = 0.1), fill = "grey70") + geom_vline(xintercept = 1, colour = "black") + theme_bw()
title <- "Distribution of AS TSS along their host genes"
for (extension in c(".png", ".pdf")) {
  ggsave(paste0(title, extension), plot = p, width = 6, height = 5, units = "in")
}

#Fig. S6B: histogram of distances between casTSS and coding TSS (2019-03-04)
df2 <- data.frame("Absolute_position" = mcols(as)$abs_pos)
title <- "Absolute distances between AS TSS and coding TSS"
p <- ggplot(df2, aes(x = Absolute_position)) + geom_histogram(bins = 100) + ggtite(title)
for (extension in c(".png", ".pdf")) {
  ggsave(paste0(title, extension), plot = p, width = 6, height = 5, units = "in")
}

### Fig. 7B (histogram of distances between casTSS and 1st/2nd exon end):
# Consider only convergent AS transcripts and their hist genes:
conv <- mcols(as)$category == "Convergent"
as_conv <- as[conv]
host_genes_conv <- host_genes[conv]
as_flipped_tss_conv <- as_flipped_tss[conv]
# Find all first and second exons:
ebg_red <- reduce(ebg_araport)
ebg_red_unl <- unlist(ebg_red)
enum_fw <- unlist(sapply(lengths(ebg_red), function(x) { seq(1, x) })) # enumerate continuous exonic intervals in both forward and reverse directions
enum_rev <- unlist(sapply(lengths(ebg_red), function(x) { seq(x, 1) }))
e_fw <- ifelse(strand(ebg_red_unl) == "+", enum_fw, enum_rev)
e_rev <- ifelse(strand(ebg_red_unl) == "+", -enum_rev, -enum_fw)
first_exons <- ebg_red_unl[as.logical(e_fw == 1 & e_rev != -1)] # first exons in genes with at least two exons
second_exons <- ebg_red_unl[as.logical(e_fw == 2)] # second exons
# Join first and second exons with host genes for convergent transcripts:
have_first_exon <- names(host_genes_conv) %in% names(first_exons)
host_genes_conv_ex <- host_genes_conv[have_first_exon] # n = 1355
as_flipped_tss_conv_ex <- as_flipped_tss_conv[have_first_exon]
df1 <- data.frame("gene_id" = mcols(host_genes_conv_ex)$gene_id)
df2 <- data.frame("gene_id" = names(first_exons), "idx_first" = 1:length(first_exons))
df3 <- data.frame("gene_id" = names(second_exons), "idx_second" = 1:length(second_exons))
df <- left_join(df1, df2, by = "gene_id", all.x = TRUE)
df <- left_join(df, df3, by = "gene_id")
first_exons_par <- first_exons[df$idx_first]
second_exons_par <- second_exons[df$idx_second]
# Plot histogram of distances between casTSS and end of the first exon:
first_exon_gap <- pgap(as_flipped_tss_conv_ex, resize(first_exons_par, 1, "end"))
abs_from_first_exon <- ifelse(poverlapsGRanges(flank(first_exon_gap, 1), as_flipped_tss_conv_ex) & width(first_exon_gap) > 0, -width(first_exon_gap), width(first_exon_gap))
ttl <- "Absolute distance between casTSS and end of the first exon"
p <- ggplot(data.frame("Dist" = abs_from_first_exon), aes(x = Dist)) + geom_histogram(bins = 100) + 
  ggtitle(ttl) + xlab("Distance from exon end, bp") + ylab("Number of transcripts")
for (ext in c("pdf", "png")) {
  ggsave(paste(ttl, ext, sep = "."), plot = p, width = 7, height = 7, units = "in")
}
# Plot histogram of distances between casTSS and end of the second exon:
second_exon_gap <- pgap(as_flipped_tss_conv_ex, resize(second_exons_par, 1, "end"))
abs_from_second_exon <- ifelse(poverlapsGRanges(flank(second_exon_gap, 1), as_flipped_tss_conv_ex) & width(second_exon_gap) > 0, -width(second_exon_gap), width(second_exon_gap))
ttl <- "Absolute distance between casTSS and end of the second exon"
p <- ggplot(data.frame("Dist" = abs_from_second_exon), aes(x = Dist)) + geom_histogram(bins = 100) + 
  ggtitle(ttl) + xlab("Distance from exon end, bp") + ylab("Number of transcripts")
for (ext in c("pdf", "png")) {
  ggsave(paste(ttl, ext, sep = "."), plot = p, width = 7, height = 7, units = "in")
}

### Fig. S6C (nucleosome count between casTTT and coding TSS:
# Load nucleosome positions (PlantDHS):
nps <- read.table("Ath_leaf_NPS.gff", sep = "\t", header = FALSE) # download from http://plantdhs.org/static/download/Ath_leaf_NPS.gff.gz
nps <- GRanges(seqnames=sub("Chr", "", nps$V1), ranges=IRanges(nps$V4, end=nps$V5), seqinfo=seqinfo(genes_araport_adj))
nps_c <- resize(nps, width = 1, fix = "center") # centers of nucleosomes
# Enumerate all nucleosomes:
hits2 <- findOverlaps(genes_araport_adj, nps_c)
genes_p <- genes_araport_adj[queryHits(hits2)]
nps_p <- nps_c[subjectHits(hits2)]
strand(nps_p) <- strand(genes_p)
rlen <- runLength(Rle(queryHits(hits2)))
enum_fw <- unlist(sapply(rlen, function(x) { seq(1, x) }))
enum_rev <- unlist(sapply(rlen, function(x) { seq(x, 1) }))
mcols(nps_p)$enum <- ifelse(strand(nps_p) == "+", enum_fw, enum_rev)
mcols(nps_p)$gene_id <- names(genes_p)
# Find the nucleosomes between coding TSS and casTSS (can be multiple hits if a nucleosome was shared between 2 or more genes):
names(as_flipped_tss_conv) <- names(host_genes_conv)
mcols(as_flipped_tss_conv) <- NULL
hits_up <- follow(as_flipped_tss_conv, nps_p, select = "all")
conv_up <- as_flipped_tss_conv[queryHits(hits_up)]
nps_up <- nps_p[subjectHits(hits_up)]
valid_up <- names(conv_up) == mcols(nps_up)$gene_id
nums_up <- mcols(nps_up)$enum[valid_up]
# Draw the figure:
ttl <- "Convergent transcripts - Upstream nucleosomes count"
pdf(paste0(ttl, ".pdf"))
plot(table(nums_up), xlab = "Upstream nucleosomes count", ylab = "Number of transcripts"); dev.off()
png(paste0(ttl, ".png"))
plot(table(nums_up), xlab = "Upstream nucleosomes count", ylab = "Number of transcripts"); dev.off()

### Fig. S6G (boxplot of FPKM values in CAS host genes vs other genes:
# Count plaNET-Seq tags on all genes:
genes_fpkm_wt <- round(getOverlappingScores(genes_araport_adj, planet_data_norm1M[1], value = "count_matrix")/ width(genes_araport_adj) * 1000, 6)
mcols(genes_araport_adj)$FPKM_plaNET_WT <- as.numeric(genes_fpkm_wt)
# Filter out silent genes (with FPKM < 1):
genes_expr <- genes_araport_adj[mcols(genes_araport_adj)$FPKM_plaNET_WT >= 1]
# Classify transcribed genes into host genes for convergent transcripts and other genes:
mcols(genes_expr)$host_conv <- mcols(genes_expr)$gene_id %in% mcols(host_genes_conv)$gene_id # 1256 TRUE, 19449 FALSE
# Draw the boxplot:
df <- data.frame("Expression" = log2(mcols(genes_expr)$FPKM_plaNET_WT), "Group" = ifelse(mcols(genes_expr)$host_conv, "With conv", "Without conv"))
p <- ggplot(df, aes(x = Group, y = Expression)) + geom_boxplot() + ylab("Log2 FPKM") + ggtitle("Expressed genes (plaNET-Seq WT FPKM >= 1)")
for (ext in c("pdf", "png")) {
  ggsave(paste("Expression levels of genes with or without convergent transcripts", ext, sep = "."), plot = p, width = 6, height = 6, units = "in")
}
# Mann-Whitney p-value:
wilcox.test(x = df[df$Group == "With conv", "Expression"], y = df[df$Group == "Without conv", "Expression"])$p.value # 5.527256e-38

### Fig. 7D (overlap of casTSS with chromatin states):
# Group ChromHMM states from PCSD database (http://systemsbiology.cau.edu.cn/chromstates/download.php) into 5 chromatin segments:
bed_dir <- "." # change to your directory with BED files downloaded from PCSD
prom_files <- paste0("At_segments_S", c(13, 15:21), ".bed") # transcription initiation
trans_prom_early_files <- paste0("At_segments_S", 22:23, ".bed") # transition from initiation to early elongation
early_files <- paste0("At_segments_S", 24:26, ".bed") # early elongation
late_files <- paste0("At_segments_S", c(3:12, 27:28), ".bed") # late elongation
pa_files <- paste0("At_segments_S", 1:2, ".bed") # transcription termination
all_files <- list(prom_files, trans_prom_early_files, early_files, late_files, pa_files)
all_data <- vector("list", length(all_files))
types <- c("Prom", "Prom2Early", "Early", "Late", "pA")
names(all_data) <- types
colors <- brewer.pal(8, "Dark2")[c(1, 6, 2, 4, 3)]

for (i in seq_along(all_data)) {
  files <- all_files[[i]]
  data <- vector("list", length(files))
  for (j in seq_along(files)) {
    bed_file <- files[[j]]
    bed <- import(file.path(bed_dir, bed_file), format = "BED")
    seqinfo(bed, new2old = c(1:5, NA, NA), pruning.mode = "coarse") <- seqinfo(genes_araport_adj)
    mcols(bed)$name <- NULL
    data[[j]] <- bed
  }
  gr <- reduce(Reduce(c, data), min.gapwidth = 100)
  mcols(gr)$name <- names(all_data)[[i]]
  mcols(gr)$itemRgb <- colors[[i]]
  all_data[[i]] <- gr
}

segments <- Reduce(c, all_data)
segments <- sortSeqlevels(segments)
segments <- sort(segments)

# Find overlaps between casTSS and chromatin segments:
mcols(as_conv)$hmm <- "No"
hits <- findOverlaps(as_flipped_tss_conv, segments)
segm_par <- segments[subjectHits(hits)]
segm_name <- mcols(segm_par)$name
mcols(as_conv)$hmm[queryHits(hits)] <- segm_name
# Calculate the percent of convergent transcripts starting in each segment:
conv_over <- table(mcols(as_conv)$hmm) / length(as_conv)
conv_over <- conv_over[c(5, 6, 1, 2, 4)]
random_over <- vector("numeric", length(types))
names(random_over) <- types
# Calculate also the overlap of control intervals (first 50% of each host gene) with the segments:
ctrl <- resize(host_genes_conv, width(host_genes_conv) * 0.5, "start")
for (i in seq_along(types)) {
  type <- types[[i]]
  gr <- segments[mcols(segments)$name == type]
  overlap <- sum(width(reduce(GenomicRanges::intersect(gr, ctrl, ignore.strand = TRUE))))
  overlap_norm <- overlap / sum(width(reduce(ctrl)))
  random_over[[i]] <- overlap_norm
}
# Draw barplot showing the observed vs expected frequency of overlaps:
df1 <- data.frame("Freq" = as.numeric(as_over), "Segments" = names(as_over), "Type" = "Observed")
df2 <- data.frame("Freq" = unname(random_over), "Segments" = names(random_over), "Type" = "Expected")
df <- rbind(df1, df2)
df$Type <- factor(df$Type, levels = c("Expected", "Observed"))
df$Segments <- factor(df$Segments, levels = types)
title <- "Fold enrichment of Convergent TSS in HMM segments"
p <- ggplot(df, aes(x = Segments, y = Freq, fill = Type)) + geom_bar(stat = "identity", position = "dodge") + ggtitle(title)
for (ext in c(".png", ".pdf")) {
  ggsave(paste0(title, ext), plot = p, width = 6, height = 5, units = "in")
}






