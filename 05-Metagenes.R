### The following code allows to reproduce all metagene plots which were included into the paper

library(rtracklayer)
library(ggplot2)
set.seed(42)

# Load adjusted Araport11 genes (see 04-Adjustment_Araport11.R for details):
genes_araport_adj <- readRDS("genes_araport_adj.RDS")

# Extract protein-coding genes (for normalization of plaNET-Seq and other tracks):
genes_npcd <- genes_araport_adj[mcols(genes_araport_adj)$tx_type == "mRNA" & seqnames(genes_araport_adj) %in% 1:5] # nuclear protein-coding genes
genes_npcd_ext <- suppressWarnings(trim(resize(genes_npcd, width(genes_npcd) + 100, "end"))) # extend by 100 bp upstream
genes_npcd_ext <- suppressWarnings(trim(resize(genes_npcd_ext, width(genes_npcd_ext) + 500, "start")))  # extend by 500 bp downstream to capture pA peaks in plaNET-Seq data
genes_npcd_m <- reduce(genes_npcd_ext) # merge overlapping intervals

# Load the original Araport11 annotation:
txdb_araport <- makeTxDbFromGFF("Araport11_GFF3_genes_transposons.201606.gff.gz")
ebg_araport <- exonsBy(txdb_araport, by = "gene") # exons grouped by gene
seqinfo(ebg_araport, new2old = as.integer(c(1:5, 7, 6))) <- seqinfo(genes_araport_adj)
exons <- unlist(ebg_araport)
ebg_npcd <- ebg_araport[names(ebg_araport) %in% names(genes_npcd)] # exons in nuclear protein-coding genes
#exons_npcd <- unlist(ebg_npcd)

# Load custom functions:
r_dir <- "D:/SCIENCE-4/My R scripts"
scripts <- c("batchReadTrackData.R", "metageneMatrix.R", "drawMetagenePlot.R", "normalizeGR.R", "removeFirstAndLastExons.R", 
             "mergeBedgraphs.R", "findFirstNucleosome.R", "findMatchedControl.R", "getOverlappingScores.R")
for (script in scripts) { source(file.path(r_dir, script)) }

# Load plaNET-Seq data (merged biological replicates):
planet_dir <- "." # change to the directory containing merged plaNET-Seq Bedgraph files obtained from 02-Postprocessing_plaNET-Seq.R
planet_files <- list.files(planet_dir, pattern = "merged_fw_rev.bedgraph.gz$")
planet_data <- batchReadTrackData(planet_files, dir = planet_dir, format = "bedGraph", seqinfo = seqinfo(genes_araport_adj))
names(planet_data) <- paste0("plaNET_", sub("_merged_fw_rev.bedgraph.gz", "", names(planet_data)))
planet_data <- endoapply(planet_data, normalizeGR, by = genes_npcd_m) # normalize track to 1M tags in nuclear protein-coding genes

# Load and normalize the remapped pNET-Seq data (merged replicates):
pnet_dir <- "." # change to the relevant directory
pnet_files <- list.files(pnet_dir, pattern = "merged_fw_rev.bedgraph.gz$")
pnet_data <- batchReadTrackData(pnet_files, dir = pnet_dir, format = "bedGraph", seqinfo = seqinfo(genes_araport_adj))
names(pnet_data) <- paste0("pNET_", sub("_merged_fw_rev.bedgraph.gz", "", names(pnet_data)))
pnet_data <- endoapply(pnet_data, normalizeGR, by = genes_npcd_m)

# Load and normalize the remapped GRO-Seq data (merged replicates; see 03-Alignment_GRO-Seq_RNA-Seq_DR-Seq):
gro_dir <- "." # change to the relevant directory
gro_files <- list("GRO_Liu2018_nrpd1e1_merged_fw_rev.bedgraph.gz", "GRO_Zhu2018_merged_fw_rev.bedgraph.gz")
gro_data <- batchReadTrackData(gro_files, dir = gro_dir, format = "bedGraph", seqinfo = seqinfo(genes_araport_adj))
names(gro_data) <- sub("_merged_fw_rev.bedgraph.gz", "", names(gro_data))
gro_data <- endoapply(gro_data, normalizeGR, by = genes_npcd_m)

# Load nucleosome occupancy data (PlantDHS):
nucl_dir <- "." # change to the relevant directory
nucl_file <- "Ath_leaf_NPS.bw" # download from http://plantdhs.org/static/download/Ath_leaf_DNase.bw
nucl_data <- import(file.path(nucl_dir, nucl_file), format = "BigWig")
seqinfo(nucl_data, new2old = c(1:5, NA, NA)) <- seqinfo(genes_araport_adj) # Ath_leaf_NPS.bw contains Chr1..Chr5 instead of 1..5,Mt,Pt
nucl_data <- list("Nucleosome_PlantDHS" = nucl_data)
nucl_data <- endoapply(nucl_data, normalizeGR, by = genes_npcd_m)

# Load nucleosome positions (PlantDHS):
nps <- read.table("Ath_leaf_NPS.gff", sep = "\t", header = FALSE) # download from http://plantdhs.org/static/download/Ath_leaf_NPS.gff.gz
nps <- GRanges(seqnames=sub("Chr", "", nps$V1), ranges=IRanges(nps$V4, end=nps$V5), seqinfo=seqinfo(genes_araport_adj))

# Load TSS-Seq tracks:
tss_dir <- "." # directory with "Plus" and "Minus" BigWig TSS-Seq files from wild type (GSE113677) and hen2-2 (GSE131733) 
tss_files_fw <- batchReadTrackData(list.files(tss_dir, pattern = "Plus.bw$"), dir = tss_dir, format = "BigWig")
strand(tss_files_fw) <- "+"
tss_files_rev <- batchReadTrackData(list.files(tss_dir, pattern = "Minus.bw$"), dir = tss_dir, format = "BigWig")
strand(tss_files_rev) <- "-"
# Merge Fw and Rev files:
tss_files <- mapply(function(x, y) { return(sort(c(x, y))) }, tss_files_fw, tss_files_rev, SIMPLIFY = FALSE)
names(tss_files) <- sub("_Plus.bw", "", names(tss_files))
# Merge Rep1 and Rep2:
tss_wt <- sort(c(tss_files[[1]], tss_files[[2]]))
tss_hen2 <- sort(c(tss_files[[3]], tss_files[[4]]))
tss_data <- list("TSS-Seq_WT" = tss_wt, "TSS-Seq_hen2-2" = tss_hen2)
tss_data <- endoapply(tss_data, normalizeGR, by = genes_npcd_m)


########## MAKE GENOMIC INTERVALS ##########

# 1) [TSS-500bp, TSS+500bp]:
tss_500 <- suppressWarnings(trim(resize(resize(genes_npcd, 0, "start"), 1000, "center")))
tss_500_filt <- tss_500[width(tss_500) == 1000 & countOverlaps(tss_500, genes_araport_adj) == 1] # skip windows which overlap with other annotated genes
# 2) [PAS-500bp, PAS+500bp]:
pas_500 <- suppressWarnings(trim(resize(resize(genes_npcd, 0, "end"), 1000, "center")))
pas_500_filt <- pas_500[width(pas_500) == 1000 & countOverlaps(pas_500, genes_araport_adj) == 1]
# 3) Internal exons (without the first and the last continuous exonic intervals):
int_exons <- removeFirstAndLastExons(ebg_npcd)[[3]]
int_exons_filt <- int_exons[countOverlaps(int_exons, genes_araport_adj) == 1 & countOverlaps(int_exons, exons) == 1]
int_exons_50_300 <- int_exons_filt[width(int_exons_filt) >= 50 & width(int_exons_filt) <= 300] # subset exons to width from 50 bp to 300 bp
int_exons_50_300 <- resize(int_exons_50_300, width(int_exons_50_300) - 10, "center") # trim exons by 5 bp each side to avoid possible edge effects
# 4) Introns:
ibg_araport <- psetdiff(unlist(range(ebg_npcd)), ebg_npcd) # intronic intervals grouped by gene
introns <- unlist(ibg_araport)
introns_filt <- introns[countOverlaps(introns, genes_araport_adj) == 1]
introns_50_300 <- introns_filt[width(introns_filt) >= 50 & width(introns_filt) <= 300] # subset introns to width from 50 bp to 300 bp
introns_50_300 <- resize(introns_50_300, width(introns_50_300) - 10, "center") # trim introns by 5 bp each side
# 5) Centers of first nucleosomes:
nps_first <- findFirstNucleosome(nps, genes_npcd) # for each gene, find the first nucleosome (within 500 bp downstream from the annotated TSS)
nps_first_500 <- suppressWarnings(trim(resize(nps_first, 1000, "center"))) # make 1 Kb windows around centers of the first nucleosomes
nps_first_500_filt <- nps_first_500[width(nps_first_500) == 1000 & countOverlaps(nps_first_500, genes_araport_adjust) == 1]

########## DRAW METAGENE PLOTS OF PLANET-SEQ AND PNET-SEQ TRACKS AT TSS, FIRST NUCLEOSOMES, PAS, INTERNAL EXONS AND INTRONS ##########
# (for Fig. 2A, 2C, 3C-D, 4D, S2A-C, S4C)

all_intervals <- list(list(tss_500_filt, "TSS 500 bp", TRUE, 200, "Window (1 Kb) centered at TSS (5 bp bins)", TRUE),
                      list(pas_500_filt, "PAS 500 bp", TRUE, 200, "Window (1 Kb) centered at PAS (5 bp bins)", TRUE),
                      list(nps_first_500_filt, "First nucl 500 bp", TRUE, 200, "Window (1 Kb) centered at first nucleosome (5 bp bins)", TRUE),
                      list(int_exons_50_300, "Exons 50bp to 300bp", TRUE, 100,  "Exons 50-300 bp scaled to 100 bins", FALSE),
                      list(introns_50_300, "Introns 50bp to 300bp", TRUE, 100, "Introns 50-300 bp scaled to 100 bins", FALSE),
                      list(introns_50_300, "Introns 50bp to 300bp unscaled", FALSE, 100, "Introns 50-300 bp unscaled (trimmed to 100 bp) anchored at start", FALSE))

# Change indexes in square brackets according to the order of samples in your planet_data and pnet_data lists:
all_data <- list(c("planet_data[1]", "PlaNET WT"),
                 c("planet_data[1:3]", "PlaNET WT Cold"),
                 c("planet_data[c(1, 4)]", "PlaNET WT Fas2"),
                 c("planet_data[c(5, 6)]", "PlaNET DMSO PlaB"),
                 c("pnet_data", "pNET Ab"))


for (i in seq_along(all_data)) {
  curr_data <- all_data[[i]]
  data <- eval(parse(text = curr_data[[1]])) # convert string to variable name
  ttl2 <- curr_data[[2]] # part 2 of the title
  message(ttl2); flush.console()
  for (j in seq_along(all_intervals)) {
    curr_int <- all_intervals[[j]]
    int <- curr_int[[1]] # current genomic windows
    ttl1 <- curr_int[[2]] # part 1 of the title
    sc <- curr_int[[3]] # scaling of genomic windows allowed (TRUE/FALSE)
    mlen <- curr_int[[4]] # ncol of the metagene matrix
    xlab <- curr_int[[5]] # label for the X axis
    vl <- curr_int[[6]] # draw vertical line at 0 (TRUE/FALSE)
    message("\t", ttl1); flush.console()
    matlist <- lapply(data, metageneMatrix, intervals = int, skip.zeros = FALSE, scaling = sc, matrix.length = mlen)
    if (isTRUE(vl)) {
      ml <- ncol(matlist[[1]])
      x.axis <- seq(-(ml/2-1), ml/2)
      vline <- 0
    } else {
      x.axis <- FALSE
      vline <- FALSE
    }
    drawMetagenePlot(matlist, x.axis = x.axis, vline = vline, title = paste0(ttl1, " (n=", nrow(matlist[[1]]), ") ", ttl2), xlabel = xlab,
                     ylim = c(0, NA), width = 8, height = 8, units = "in")
  }
}


########## DRAW METAGENE PLOTS OF SHORT VS LONG INTRONS ##########

# 1) Nucleosome occupancy of introns stratified by length (Fig. S4E):
introns_min25 <- introns[width(introns) >= 25] 
breaks = c(60, 150, 250, 500, 1000)
groups <- cut(width(introns_min25), breaks = breaks)
mcols(introns_min25)$groups <- groups
introns_trimmed <- resize(introns_min25, width(introns_min25) - 10, "center")
introns_filt <- introns_trimmed[!is.na(mcols(introns_trimmed)$groups)]

lvl <- levels(mcols(introns_filt)$groups)
ml <- vector("list", length(lvl))

for (i in seq_along(lvl)) {
  curr_lvl <- lvl[[i]]
  curr_introns <- introns_filt[mcols(introns_filt)$groups == curr_lvl]
  ml[[i]] <- metageneMatrix(signal = nucl_data[[1]], intervals = curr_introns, scaling = TRUE, matrix.length = 300, skip.zeros = FALSE)
}
names(ml) <- paste0(lvl, " (n=", lapply(ml, nrow), ")")
drawMetagenePlot(ml, title = "Introns grouped by width - Nucleosome occupancy", xlabel = "Introns scaled to 300 bins", 
                 ylabel = "Nucleosome occupancy", width = 8, height = 8, units = "in")

# 2) plaNET-Seq and pNET-Seq profiles of the first 200 bp of short vs long introns (Fig. 4F-G):
introns_60_250 <- introns_filt[mcols(introns_filt)$groups %in% lvl[2:3]] # "(60,150]" and "(150,250]"
introns_250_1000 <- introns_filt[mcols(introns_filt)$groups %in% lvl[4:5]] # "(250,500]" and "(500,1e+03]"

matlist <- lapply(planet_data[c(5, 6)], metageneMatrix, intervals = introns_60_250, skip.zeros = FALSE, scaling = FALSE, matrix.length = 200)
drawMetagenePlot(matlist, title = paste0("Short introns unscaled (first 200bp) (n=", nrow(matlist[[1]]), ") PlaNET DMSO PlaB"), 
                 xlabel = "Introns 60-250bp unscaled (trimmed to 200bp) anchored at start", ylim = c(0, NA), width = 8, height = 8, units = "in")
matlist <- lapply(planet_data[c(5, 6)], metageneMatrix, intervals = introns_250_1000, skip.zeros = FALSE, scaling = FALSE, matrix.length = 200)
drawMetagenePlot(matlist, title = paste0("Long introns unscaled (first 200bp) (n=", nrow(matlist[[1]]), ") PlaNET DMSO PlaB"), 
                 xlabel = "Introns 250-1000bp unscaled (trimmed to 200bp) anchored at start", ylim = c(0, NA), width = 8, height = 8, units = "in")
matlist <- lapply(pnet_data, metageneMatrix, intervals = introns_60_250, skip.zeros = FALSE, scaling = FALSE, matrix.length = 200)
drawMetagenePlot(matlist, title = paste0("Short introns unscaled (first 200bp) (n=", nrow(matlist[[1]]), ") pNET Ab"), 
                 xlabel = "Introns 60-250bp unscaled (trimmed to 200bp) anchored at start", ylim = c(0, NA), width = 8, height = 8, units = "in")
matlist <- lapply(pnet_data, metageneMatrix, intervals = introns_250_1000, skip.zeros = FALSE, scaling = FALSE, matrix.length = 200)
drawMetagenePlot(matlist, title = paste0("Long introns unscaled (first 200bp) (n=", nrow(matlist[[1]]), ") pNET Ab"), 
                 xlabel = "Introns 250-1000bp unscaled (trimmed to 200bp) anchored at start", ylim = c(0, NA), width = 8, height = 8, units = "in")

# 3) PlaNET-Seq profile of matched pairs of long and short introns (Fig. S4F):
# Find pairs of long-short introns matched by transcription (for each long intron, pick up a short intron from the same gene):
introns_500_1000 <- introns_filt[mcols(introns_filt)$groups == lvl[[5]]]
introns_60_500 <- introns_filt[mcols(introns_filt)$groups %in% lvl[2:4]] # introns 60-500 bp are considered "short"
dupl <- duplicated(names(introns_500_1000))
introns_500_1000_dedup <- introns_500_1000[!dupl] # to avoid repeatedly picking up the same short intron from a gene having two long introns
introns_60_500_grl <- split(introns_60_500, names(introns_60_500)) # group short introns by gene name
df1 <- data.frame("Gene_names" = names(introns_500_1000_dedup), stringsAsFactors = FALSE)
df2 <- data.frame("Gene_names" = names(introns_60_500_grl), "Index" = 1:length(introns_60_500_grl), stringsAsFactors = FALSE)
df <- left_join(df1, df2, by = "Gene_names")
not_na <- !is.na(df$Index)
introns_500_1000_dedup2 <- introns_500_1000_dedup[not_na]
introns_60_500_grl_par <- introns_60_500_grl[df$Index[not_na]] # make GRangesList of short introns parallel to long introns
random_choice <- lapply(lengths(introns_60_500_grl_par), function(x) { out <- rep(FALSE, x); idx <- sample(seq(1, x), size = 1); out[idx] <- TRUE; return(out) })
matched_introns_60_500 <- unlist(introns_60_500_grl_par)[unlist(random_choice)] # randomly choose one short introns for each long intron
# Plot plaNET WT track over unscaled long introns and matched short introns:
mat1 <- metageneMatrix(signal = planet_data[[1]], intervals = introns_500_1000_dedup2, scaling = FALSE, matrix.length = 800, 
                       skip.zeros = FALSE, skip.top.obs = TRUE, shrink = TRUE)
mat2 <- metageneMatrix(signal = planet_data[[1]], intervals = matched_introns_60_500, scaling = FALSE, matrix.length = 300, 
                       skip.zeros = FALSE, skip.top.obs = TRUE, shrink = TRUE)
mat2_ext <- cbind(mat2, matrix(nrow = nrow(mat2), ncol = (800 - 300) / 5))
ml <- list("Long introns (500-1000 bp)" = mat1, "Matched short introns (60-500 bp)" = mat2_ext)
drawMetagenePlot(ml, title = "Unscaled matched introns grouped by width - plaNET WT", xlabel = "Unscaled introns (5 bp bins) anchored as 5pSS", 
                 ylabel = "PlaNET-Seq WT", width = 10, height = 7, units = "in")


########## DRAW METAGENE PLOTS OF DIVERGENT PROMOTERS ##########

# Load coordinates of novel transcripts which were called by groHMM from plaNET-Seq data (see 06-groHMM_pipeline.R)
grohmm_novel <- readRDS("groHMM_merged_novel.RDS")
# Extract divergent transcripts (start within 500 bp from a protein-coding gene):
dnc <- subset(grohmm_novel, DNC & !Subint) # n = 917
# Find their mate genes:
genes_flipped <- genes_araport_adj
strand(genes_flipped) <- ifelse(strand(genes_flipped) == "+", "-", "+")
idx <- follow(dnc, genes_flipped)
mate_genes <- genes_araport_adj[idx] 
# Subset to nuclear protein-coding mate genes:
npcd <- mcols(mate_genes)$tx_type == "mRNA" & seqnames(mate_genes) %in% 1:5
dnc <- dnc[npcd]
mate_genes <- mate_genes[npcd]
# Deduplicate on mate genes:
dupl <- duplicated(names(mate_genes))
dnc <- dnc[!dupl]
mate_genes <- mate_genes[!dupl]

# 1) Nucleosome occupancy in the gap between DNC TSS and mate gene TSS (Fig. 6D):
between <- pgap(mate_genes, dnc, ignore.strand = TRUE) # gap has the same strandness as mate genes
win_up <- flank(between, 250)
win_down <- flank(between, 250, start = FALSE)
m1 <- metageneMatrix(signal = nucl_data[[1]], intervals = win_up, scaling = TRUE, matrix.length = 50, skip.zeros = FALSE, skip.outliers = FALSE)
m2 <- metageneMatrix(signal = nucl_data[[1]], intervals = between, scaling = TRUE, matrix.length = 50, skip.zeros = FALSE, skip.outliers = FALSE)
m3 <- metageneMatrix(signal = nucl_data[[1]], intervals = win_down, scaling = TRUE, matrix.length = 50, skip.zeros = FALSE, skip.outliers = FALSE)
matlist <- list("MNase-Seq_PlantDHS" = cbind(m1, m2, m3))
drawMetagenePlot(matlist, x.axis = seq(-49, 100), title = "Nucleosome occupancy between DNC TSS and mate gene TSS", 
                 xlabel = "Gaps scaled to 50 bins with 250 bp flanks (5 bp bins)", vline = c(0, 50), width = 8, height = 8, units = "in")

# 2) Nucleosome occupancy of DNC promoters vs matched non-DNC promoters (Fig. 6C):
# Calculate FPKM transcription values of: i) mate genes; ii) all nuclear protein-coding genes:
mate_genes_fpkm <- as.numeric(getOverlappingScores(mate_genes, planet_data[1], value = "count_matrix")) / width(mate_genes) * 1000
genes_npcd_fpkm <- as.numeric(getOverlappingScores(genes_npcd, planet_data[1], value = "count_matrix")) / width(genes_npcd) * 1000
# Find control set of nuclear protein-coding genes which have the same distribution of expression values as the mate genes: 
matched <- genes_npcd[findMatchedControl(mate_genes_fpkm[mate_genes_fpkm > 0], genes_npcd_fpkm, log = TRUE)]
# Make windows [-1.0Kb, +0.5Kb] around TSS of mate genes and control genes:
win_mate <- resize(resize(granges(mate_genes), 500, "start"), 1500, "end")
win_matched <- resize(resize(granges(matched), 500, "start"), 1500, "end")
mat1 <- metageneMatrix(signal = nucl_data[[1]], intervals = win_mate, scaling = TRUE, matrix.length = 300, skip.zeros = FALSE)
mat2 <- metageneMatrix(signal = nucl_data[[1]], intervals = win_matched, scaling = TRUE, matrix.length = 300, skip.zeros = FALSE)
matlist <- list("Mate genes" = mat1, "Matched non-DNC genes" = mat2)
drawMetagenePlot(matlist, x.axis = seq(-199, 100), title = "Nucleosome occupancy at gene TSS", 
                 xlabel = "1Kb upstream + 0.5Kb downstream from gene TSS (5 bp bins)", vline = 0, width = 8, height = 8, units = "in")

# 3) Draw pNET-Seq, GRO-Seq, TSS-Seq (WT and hen2-2) and plaNET-Seq (WT vs fas2-4) profiles around divTSS (Fig. S5A-D):
windows <- resize(resize(dnc, 0, "start"), 1000, "center")
win_flipped <- granges(windows)
strand(win_flipped) <- ifelse(strand(win_flipped) == "+", "-", "+")

cases <- list("Nucleosome occupancy" = list("nucl_data", FALSE),
              "GRO-Seq" = list("gro_data", TRUE),
              "TSS-Seq WT hen2" = list("tss_data", TRUE),
              "PlaNET-Seq WT Fas2" = list("planet_data[c(1, 4)]", TRUE),
              "pNET-Seq" = list("pnet_data", TRUE))

xlab = "1Kb windows centered at DNC TSS (5 bp bins)"

for (i in seq_along(cases)) {
  curr_name <- names(cases)[[i]]
  curr_data <- eval(parse(text = cases[[i]][[1]]))
  curr_as <- cases[[i]][[2]] # antisenseMode parameter for metageneMatrix() function (TRUE/FALSE) 
  matlist_s <- lapply(curr_data, metageneMatrix, intervals = win_flipped, scaling = TRUE, matrix.length = 200, skip.zeros = FALSE)
  drawMetagenePlot(matlist_s, x.axis = seq(-99, 100), vline = 0, title = paste0("DNC - ", curr_name, " (Fw)"), 
                   xlabel = xlab, width = 8, height = 8, units = "in")
  if (isTRUE(curr_as)) {
    matlist_as <- lapply(curr_data, metageneMatrix, intervals = win_flipped, scaling = TRUE, matrix.length = 200, skip.zeros = FALSE, antisenseMode = TRUE)
    matlist_as <- lapply(matlist_as, function(mat) { return(-mat) })
    drawMetagenePlot(matlist_as, x.axis = seq(-99, 100), vline = 0, title = paste0("DNC - ", curr_name, " (Rev)"), 
                     xlabel = xlab, width = 8, height = 8, units = "in", hline = 0)
  }
}


########## DRAW METAGENE PLOTS OF CONVERGENT PROMOTERS ##########
# Nucleosome occupancy at casTSS (Fig. 7C);
# PlaNET-Seq signal (WT vs Cold 3h vs Cold 12h) at casTSS on both sense and antisense strands (Fig. 7E);
# pNET-Seq signal at casTSS on the antisense strand (Fig. S6A);
# TSS-Seq signal (WT vs hen2-2) at casTSS on the antisense strand (Fig. S6F);

# Extract all AS transcripts (either start within the host gene, or at least overlap its 3' end):
as <- subset(grohmm_novel, over_AS & N_over_AS == 1 & !Subint)
# Flip AS transcript TSS to the coding strand:
as_flipped <- as
strand(as_flipped) <- ifelse(strand(as) == "+", "-", "+")
as_flipped_tss <- resize(as_flipped, width = 1, fix = "end")
# Find host genes:
host_genes <- genes_araport_adj[findOverlaps(as_flipped, genes_araport_adj, select = "first")]
# Relative positions of AS transcripts within their host genes:
abs_pos <- width(pgap(resize(host_genes, 1, "start"), as_flipped_tss))
rel_pos <- abs_pos / width(host_genes) # can exceed 1.0
# Subset to convergent transcripts (start within the first 50% of the host gene width):
convergent <- rel_pos <= 0.5
cas_tss <- as_flipped_tss[convergent] # n = 1497 (flipped to the strand of the host gene)
# Make 1Kb windows around casTSS:
windows_fw <- granges(resize(cas_tss, width = 1000, fix = "center")) # strand of the host gene
windows_rev <- windows_fw
strand(windows_rev) <- ifelse(strand(windows_fw) == "+", "-", "+") # strand of the AS transcript
# Filter windows for no overlap with other known genes on both strands:
no_overlap <- countOverlaps(windows_fw, genes_araport_adj) == 1 & countOverlaps(windows_rev, genes_araport_adj) == 0
windows_filt <- windows_fw[no_overlap] # n = 1404
# Draw metagene plots:
xlab = "1Kb windows centered at convergent transcript start"

cases <- list("Nucleosome density" = list("nucl_data", FALSE),
              "TSS-Seq WT hen2" = list("tss_data", TRUE),
              "PlaNET-Seq WT Cold" = list("planet_data[1:3]", TRUE),
              "pNET-Seq" = list("pnet_data", TRUE))

for (i in seq_along(cases)) {
  curr_name <- names(cases)[[i]]
  curr_data <- eval(parse(text = cases[[i]][[1]]))
  curr_as <- cases[[i]][[2]] # antisenseMode parameter for metageneMatrix() function (TRUE/FALSE) 
  matlist_s <- lapply(curr_data, metageneMatrix, intervals = windows_filt, scaling = TRUE, matrix.length = 200, skip.zeros = FALSE)
  drawMetagenePlot(matlist_s, x.axis = seq(-99, 100), vline = 0, title = paste0("Convergent transcripts - ", curr_name, " (Fw)"), 
                   xlabel = xlab, width = 8, height = 8, units = "in")
  if (isTRUE(curr_as)) {
    matlist_as <- lapply(curr_data, metageneMatrix, intervals = windows_filt, scaling = TRUE, matrix.length = 200, skip.zeros = FALSE, antisenseMode = TRUE)
    matlist_as <- lapply(matlist_as, function(mat) { return(-mat) })
    drawMetagenePlot(matlist_as, x.axis = seq(-99, 100), vline = 0, title = paste0("Convergent transcripts - ", curr_name, " (Rev)"), 
                     xlabel = xlab, width = 8, height = 8, units = "in", hline = 0)
  }
}


########## THE DISTRIBUTION OF CHROMATIN SEGMENTS ALONG PROTEIN-CODING GENES (FIG. S6E) ##########

# Load chromatin states from PCSD database combined into 5 biologically relevant gene segments:
segments <- readRDS("PCSD_gene_segments.RDS") # this file was generated in the 06-groHMM_pipeline.R

# Find nuclear protein-coding genes with no overlap with any other known gene with 500 bp margins on each side:
genes_npcd_ext_1kb <- suppressWarnings(trim(resize(genes_npcd, width(genes_npcd) + 1000, "center")))
no_over <- countOverlaps(genes_npcd_ext_1Kb, genes_araport_adj) == 1
# Choose nuclear protein-coding expressed above 1 FPKM (based on plaNET-Seq WT data) and having length 1-5 Kb:
genes_good <- genes3_npcd[no_over & genes_npcd_fpkm >= 1 & width(genes_npcd) >= 1000 & width(genes_npcd) <= 5000] # n = 10291

# Make "0-1" tracks from HMM segments (the presence of a segment is denoted by interval with score 1):
segm_names <- c("Prom", "Prom2Early", "Early", "Late", "pA")
segm_tracks <- vector("list", length(segm_names))

for (i in seq_along(segm_tracks)) {
  segm_name <- segm_names[[i]]
  gr <- segments[mcols(segments)$name == segm_name]
  score(gr) <- 1
  segm_tracks[[i]] <- gr
}

# Compute matrices:
a <- lapply(segm_tracks, metageneMatrix, intervals = flank(genes_good, 250), scaling = TRUE, matrix.length = 50, skip.zeros = FALSE, skip.outliers = FALSE)
b <- lapply(segm_tracks, metageneMatrix, intervals = genes_good, scaling = TRUE, matrix.length = 300, skip.zeros = FALSE, skip.outliers = FALSE)
c <- lapply(segm_tracks, metageneMatrix, intervals = flank(genes_good, 250, start = FALSE), scaling = TRUE, matrix.length = 50, skip.zeros = FALSE, skip.outliers = FALSE)
# Reorder:
matlist <- vector("list", length(segm_names))
names(matlist) <- segm_names
for (i in seq_along(matlist)) {
  matlist[[i]] <- cbind(a[[i]], b[[i]], c[[i]])
}

# Draw the plot:
drawMetagenePlot(matlist, title = "PCSD segments in nuclear protein-coding genes FPKM1 1-5Kb (n = 10291)", x.axis = seq(-49, 350), vline = c(0, 300), 
                 xlabel = "[TSS, TTS] scaled to 300 bins with 250 bp flanks", width = 8, height = 8, units = "in")
