# This script adjusts borders of Araport11 genes using data from TIF-Seq, TSS-Seq and DR-Seq experiments;
# TIF-Seq, or Transcript Isoform Sequencing, was originally developed by Pelechano et al., 2014 (PMID 24967623) in budding yeast;
# This method allows for simultaneous mapping of 5' and 3'-ends of the same RNA molecule (capped and polyadenylated);
# I.e first bases of read 1 and read 2 in paired-end sequencing denote borders of the transcript isoform (both TSS and PAS);
# Recently, we for the first time applied TIF-Seq to Arabidopsis Col-0 seedlings (manuscript co-submitted);
# TSS-Seq is a simplification of TIF-Seq protocol which resolves only 5'-ends of polyadenylated transcripts. In effect, it is similar to CAGE;
# Our TSS-Seq dataset for Col-0 seedlings was published in Nielsen et al., 2019 (PMID 30707695) and can be downloaded from GSE113677;
# DR-Seq, or Direct RNA sequencing, allows to resolve PAS sites;
# We combined DR-Seq datasets from Schurch et al., 2014 (PMID 24722185) and Sherstnev et al., 2012 (PMID 22820990), as detailed in 03-Alignment_GRO-Seq_RNA-Seq_DR-Seq.sh;


min_tif_score <- 2 # use TIF-Seq clusters with at least 2 supporting read pairs
min_tss_score <- 5 # use TSS-Seq clusters with the score at least 5
min_pas_score <- 5 # the same for DR-Seq clusters
min_tif_gene_overlap <- 0.8 # a valid TIF-Seq cluster should overlap Araport11 gene by at least 80% of its length
min_abs_tif_cov <- 5 # threshold values for absolute and relative TIF-Seq coverage are used to adjust Araport11 genes by TIF-Seq profile alone
min_rel_tif_cov <- 0.1 # (i.e. when TSS-Seq and/or DR-Seq clusters are occassionally missing for given gene)
max_extension <- 500 # extension of gene borders is limited to 0.5 Kb
ncrna_max_shrink <- 50 # non-coding RNAs are allowed to shrink by no more than 50 bp...
long_ncrna <- 300 # ... and only given that their annotated length was not less than 300 bp
# (protein-coding genes are allowed to shrink along their 5' and 3' UTRs; they cannot shrink beyond the borders of the longest annotated CDS)

library(SummarizedExperiment)
library(rtracklayer)
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb_tair <- TxDb.Athaliana.BioMart.plantsmart28 # Ensembl Plants 28 (very similar to TAIR10)

# Araport11 annotation file was downloaded from https://www.arabidopsis.org/download_files/Genes/Araport11_genome_release/Araport11_GFF3_genes_transposons.201606.gff.gz
txdb_araport <- makeTxDbFromGFF("Araport11_GFF3_genes_transposons.201606.gff.gz")
genes <- genes(txdb_araport, columns = c("gene_id", "tx_type")) # n = 38194
seqinfo(genes, new2old = as.integer(c(1:5, 7, 6))) <- seqinfo(txdb_tair) # fix chromosome names
mcols(genes)$tx_type <- unname(unlist(mcols(genes)$tx_type))
mcols(genes)$thick <- granges(genes) # save the original Araport11 coordinates in mcols

# Add dummy "cds" to mcols (determined how much the gene can shrink during the adjustment procedure):
genes_shr <- genes
long <- width(genes) >= long_ncrna
genes_shr[long] <- resize(genes[long], width = width(genes[long]) - ncrna_max_shrink * 2, fix = "center")
mcols(genes)$cds <- granges(genes_shr)
# Change dummy "cds" to real CDS for protein-coding genes:
cds_grl <- cdsBy(txdb_araport, by = "gene")
cds_gr <- unlist(range(cds_grl))
mcols(genes)$cds[names(genes) %in% names(cds_gr)] <- granges(cds_gr)
# Compute maximal extension borders:
max_ext <- granges(suppressWarnings(trim(resize(genes, width = width(genes) + max_extension * 2, fix = "center"))))
# Add upstream and downstream intervals where the adjustment of respective gene borders is allowed:
mcols(genes)$up <- pgap(resize(max_ext, width = 0, fix = "start"), mcols(genes)$cds)
mcols(genes)$down <- pgap(mcols(genes)$cds, resize(max_ext, width = 0, fix = "end"))

# Load TIF-Seq data:
tif <- import("TIF-Seq_clusters_WT.bed", format = "BED", seqinfo = seqinfo(txdb_tair))
mcols(tif)$score <- as.integer(mcols(tif)$name)
mcols(tif)$name <- NULL
tif <- tif[mcols(tif)$score >= min_tif_score]

# Load TSS-Seq data (see Nielsen et al., 2019 - PMID 30707695 for details on TSS-Seq data processing and cluster calling by CAGEfightR):
tss <- rowRanges(readRDS("TSS-Seq_Col-0.RDS"))
tss <- tss[mcols(tss)$score >= min_tss_score]

# Load PAS data from DR-Seq (DR-Seq clusters were called by CAGEfightR in the same way as TSS-Seq data; see Nielsen et al., 2019 for details):
pas <- rowRanges(readRDS("DR-Seq_Col-0_cov2.RDS"))
pas <- pas[mcols(pas)$score >= min_pas_score]

# Skip intergenic TIF-Seq clusters which do not overlap any gene:
tif <- tif[tif %over% genes]

# Find overlaps between TIF-Seq clusters and genes:
hits <- findOverlaps(tif, genes)
tif_par <- tif[queryHits(hits)]
genes_par <- genes[subjectHits(hits)]

# Find valid overlaps (at least 80% of the gene covered):
valid <- width(pintersect(genes_par, tif_par)) / width(genes_par) >= min_tif_gene_overlap

# Detect ambiguous clusters which have valid overlaps with more than one gene: 
one_valid <- as.logical(tapply(valid, queryHits(hits), function(x) { sum(x) == 1 }))
one_valid_ext <- rep(one_valid, times = runLength(Rle(queryHits(hits)))) # Extend 'one_valid' to the length of 'valid'
final_valid <- valid & one_valid_ext

# Subset the overlaps to 'final_valid': 
tif_par2 <- tif_par[final_valid]
genes_par2 <- genes_par[final_valid]
hits2 <- hits[final_valid]

# Combine TIF-Seq clusters with their relevant genes:
tif_spl <- split(tif_par2, subjectHits(hits2)) # group valid clusters by the gene
genes_spl <- split(genes_par2, subjectHits(hits2))

out <- genes_spl
pb <- txtProgressBar(min = 1, max = length(out), style = 3)
for (i in seq_along(tif_spl)) {
  setTxtProgressBar(pb, i)
  curr_tif <- tif_spl[[i]]
  curr_coord <- genes_spl[[i]][1]
  mc <- mcols(curr_coord) # save mcols (otherwise it is lost during coordinate adjustment)
  # First try to extend the gene by TIF-Seq coverage alone:
  cov <- bindAsGRanges(score = coverage(curr_tif, weight = "score"))
  cov <- cov[score(cov) > 0] # skip zero coverage outside of the current gene
  above <- score(cov) >= min_abs_tif_cov & score(cov) / max(score(cov)) >= min_rel_tif_cov
  if (any(above)) {
    tif_range <- range(cov[above])
    strand(tif_range) <- strand(curr_coord)
    start_tif <- resize(tif_range, 0, "start")
    end_tif <- resize(tif_range, 0, "end")
    if (start_tif %within% mc$up) {
      curr_coord <- punion(start_tif, resize(curr_coord, 0, "end"), fill.gap = TRUE)
    }
    if (end_tif %within% mc$down) {
      curr_coord <- punion(resize(curr_coord, 0, "start"), end_tif, fill.gap = TRUE)
    }
  }
  # Then try to adjust gene starts by available TSS clusters from TSS-Seq:
  rel_tss <- tss[tss %over% resize(curr_tif, 1, "start")] # find TSS which overlap with relevant TIF-Seq clusters
  tss_summits <- rel_tss
  ranges(tss_summits) <- mcols(tss_summits)$thick # Extract TSS summits
  tss_summits <- tss_summits[tss_summits %within% mc$up] # filter TSS summits by the allowed upstream interval
  if (length(tss_summits) > 1) {
    max_tss_score <- max(mcols(tss_summits)$score)
    tss_summits <- tss_summits[mcols(tss_summits)$score == max_tss_score][1] # choose TSS with the highest score
  }
  if (length(tss_summits) == 1) { # if there is a valid TSS summit
    curr_coord <- punion(tss_summits, resize(curr_coord, 0, "end"), fill.gap = TRUE)
  }
  # Try to adjust gene ends by available PAS clusters from DR-Seq:
  rel_pas <- pas[pas %over% resize(curr_tif, 1, "end")]
  pas_summits <- rel_pas
  ranges(pas_summits) <- mcols(pas_summits)$thick
  pas_summits <- pas_summits[pas_summits %within% mc$down]
  if (length(pas_summits) > 1) {
    max_pas_score <- max(mcols(pas_summits)$score)
    pas_summits <- pas_summits[mcols(pas_summits)$score == max_pas_score][1]
  }
  if (length(pas_summits) == 1) {
    curr_coord <- punion(resize(curr_coord, 0, "start"), pas_summits, fill.gap = TRUE)
  }
  mcols(curr_coord) <- mc # restore mcols
  out[[i]] <- curr_coord
}

out <- unlist(out)
names(out) <- mcols(out)$gene_id
mcols(out)$source <- "Not_adj"
changed_coord <- start(out) != start(mcols(out)$thick) | end(out) != end(mcols(out)$thick)
sum(changed_coord) # n = 10464
mcols(out)$source[changed_coord] <- "Adjusted"

# Combine with genes which were not covered by valid TIF-Seq clusters:
not_adj <- genes[!(mcols(genes)$gene_id %in% mcols(out)$gene_id)]
mcols(not_adj)$source <- "Araport11"
genes_araport_adj <- sort(c(out, not_adj))

# Save the results:
saveRDS(genes_araport_adj, "genes_araport_adj.RDS")
export(genes_araport_adj, "genes_araport_adj.bed", format = "BED")



