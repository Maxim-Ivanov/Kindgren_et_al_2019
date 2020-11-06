# This script was used to call TSS clusters from TSS-Seq data and PAS clusters from Direct RNA-Seq (DR-Seq) data;
# The resultant coordinates of TSS and PAS clusters are used in 04-Adjustment_Araport11.R and 08-Readthrough_distance.R scripts;

library(CAGEfightR) # version 1.2.0
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28

r_dir <- "." # change to the folder where you saved the custom functions from https://github.com/Maxim-Ivanov/Kindgren_et_al_2019
scripts <- c("batchReadTrackData.R", "expandGRtoUnitWidth.R")
for (script in scripts) { source(file.path(r_dir, script)) }

##### PART 1: Call TSS clusters (i.e. promoters) from TSS-Seq data on wild type and hen2-2 samples #####

### 1) Call TSS clusters in wild type Col-0 (samples GSM3369630	and GSM3369631 from study GSE119304; Kindgren at al., 2018 - PMID 30385760):
# Load TSS-Seq hen2-2 data:
tss_dir <- "." # change to the directory where you downloaded BigWig files for samples GSM3369630 and GSM3369631
tss_data_fw <- batchReadTrackData(list.files(tss_dir, pattern = "Plus.bw$"), dir = tss_dir, format = "BigWig")
strand(tss_data_fw) <- "+"
tss_data_rev <- batchReadTrackData(list.files(tss_dir, pattern = "Minus.bw$"), dir = tss_dir, format = "BigWig")
strand(tss_data_rev) <- "-"
# Merge Fw and Rev files:
tss_data <- mapply(function(x, y) { return(sort(c(x, y))) }, tss_data_fw, tss_rev, SIMPLIFY = FALSE)
names(tss_data) <- c("wt_rep1", "wt_rep2")
# Filter by minimal coverage:
tss_data <- lapply(tss_data, function(gr) { return(gr[score(gr) >= 3]) })
# Expand to unit width and save as BigWig (input requirements of CAGEfightR):
tss_data_uw <- lapply(tss_data, expandGRtoUnitWidth)

for (i in seq_along(tss_data_uw)) {
  data <- tss_data_uw[[i]]
  name <- names(tss_data_uw)[[i]]
  export(data[strand(data) == "+"], paste0(name, "_unitWidth_Plus.bw"), format = "BigWig")
  export(data[strand(data) == "-"], paste0(name, "_unitWidth_Minus.bw"), format = "BigWig")
}

# Call TSS clusters by CAGEfightR:
bw_plus_filenames <- list.files(".", pattern = "unitWidth_Plus.bw$")
bw_minus_filenames <- list.files(".", pattern = "unitWidth_Minus.bw$")
bw_plus <- BigWigFileList(bw_plus_filenames)
bw_minus <- BigWigFileList(bw_minus_filenames)
sample_names <- sub("TSS-Seq_", "", sub("_unitWidth_Plus.bw", "", bw_plus_filenames))
names(bw_plus) <- sample_names
names(bw_minus) <- sample_names
design_matrix <- data.frame("Name" = sample_names, "BigWigPlus" = bw_plus_filenames, "BigWigMinus" = bw_minus_filenames, "Genotype" = "wt", row.names = sample_names)
ctss <- quantifyCTSSs(plusStrand = bw_plus, minusStrand = bw_minus, design = design_matrix, genome = seqinfo(genes_araport_adj))
tss <- quickTSSs(ctss)
tss <- rowRanges(subsetBySupport(tss)) # subset to clusters with have nonzero number of tags in both replicates
# Save the results:
saveRDS(tss, "TSS_clusters_WT.RDS") # (this file is used in 04-Adjustment_Araport11.R script)

### 1) Call TSS clusters in Col-0 containing hen2-2 mutation ():
# Repeat the code above using samples GSM3814858 and GSM3814859 from study GSE131733;
# Save the results as "TSS_clusters_hen2-2.RDS" (this file is used in 08-Readthrough_distance.R script)


##### PART 2: Call PAS clusters (i.e. terminators) from DR-Seq data #####
# The Direct RNA-Seq data are available from Sherstnev et al., 2012 (PMID 22820990) and Schurch et al. 2014 (PMID 24722185)
# For their remapping, see 03-Alignment_GRO-Seq_RNA-Seq_DR-Seq.sh script;

# Load Bedgraph files:
bg_dir <- "." # change to the directory containing DR-Seq Bedgraph files returned by 03-Alignment_GRO-Seq_RNA-Seq_DR-Seq.sh
bg_filenames <- c(list.files(bg_dir, pattern = "DRseq_Schurch2014_rep._fw_rev\\.bedgraph\\.gz$"), "DRseq_Sherstnev2012_merged_fw_rev.bedgraph.gz")
pas_data <- batchReadTrackData(file.path(bg_dir, bg_files), format = "bedGraph", seqinfo = seqinfo(txdb))
# Filter by minimal coverage:
pas_data <- lapply(pas_data, function(gr) { return(gr[score(gr) >= 2]) })
# Expand to unit width and save as BigWig:
pas_data_uw <- lapply(pas_data, expandGRtoUnitWidth)

for (i in seq_along(pas_data_uw)) {
  data <- pas_data_uw[[i]]
  name <- names(pas_data_uw)[[i]]
  export(data[strand(data) == "+"], paste0(name, "_unitWidth_Plus.bw"), format = "BigWig")
  export(data[strand(data) == "-"], paste0(name, "_unitWidth_Minus.bw"), format = "BigWig")
}

# Load BigWig files to CAGEfightR:
bw_plus_filenames <- list.files(".", pattern = "fw_ctss_cov2.bw")
bw_minus_filenames <- list.files(".", pattern = "rev_ctss_cov2.bw")
bw_plus <- BigWigFileList(bw_plus_filenames)
bw_minus <- BigWigFileList(bw_minus_filenames)
sample_names <- sub("DR-Seq_", "", sub("_unitWidth_Plus.bw", "", bw_plus_filenames))
names(bw_plus) <- sample_names
names(bw_minus) <- sample_names

# Call PAS clusters by CAGEfightR:
design_matrix <- data.frame("Name" = sample_names, "BigWigPlus" = bw_plus_filenames, "BigWigMinus" = bw_minus_filenames, "Genotype" = "wt", row.names = sample_names)
cpas <- quantifyCTSSs(plusStrand = bw_plus, minusStrand = bw_minus, design = design_matrix, genome = seqinfo(txdb), tileWidth = sum(seqlengths(txdb))) # quantify all tag clusters (TCs)
pas <- quickTSSs(cpas) # call candidate PAS
pas <- subsetBySupport(pas) # remove TCs detected in only one replicate
saveRDS(pas, "PAS_clusters_WT.RDS") # this file is used in both 04-Adjustment_Araport11.R and 08-Readthrough_distance.R scripts

