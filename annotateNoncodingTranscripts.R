# tx = GRanges object containing coordinates of transcripts to annotate by overlap with the gene model;
# genes = GRanges object containing the gene model (coordinates of all known genes);
# The function returns a logical matrix;

annotateNoncodingTranscripts  <- function(tx, genes) {
  library(GenomicRanges)
  tx <- granges(tx)
  genes <- granges(genes)
  if (!identical(seqinfo(tx), seqinfo(genes))) { stop("Ensure identical seqinfo!") }
  # Prepare the output data frame:
  mat <- matrix(data = FALSE, nrow = length(tx), ncol = 8, dimnames = list(NULL, c("Genic", "Up", "Down", "Divergent", "Convergent", "TTS_AS", "Distal_AS", "Intergenic")))
  # Find strong overlaps on the sense strand (either transcript or gene is covered by at least 50% width):
  over_s <- tx %over% genes
  hits <- findOverlaps(tx, genes)
  tx_par <- tx[queryHits(hits)]
  genes_par <- genes[subjectHits(hits)]
  overlap <- width(pintersect(tx_par, genes_par))
  over_tx <- overlap / width(tx_par)
  over_genes <- overlap / width(genes_par)
  strong_overlap <- over_tx >= 0.5 | over_genes >= 0.5
  mat[over_s, "Genic"] <- as.logical(unlist(tapply(strong_overlap, queryHits(hits), any)))
  # Detect weak upstream and downstream overlaps on the same strand:
  weak_overlap <- !strong_overlap
  a <- start(tx_par) <= start(genes_par) & end(tx_par) < end(genes_par)
  b <- start(tx_par) > start(genes_par) & end(tx_par) >= end(genes_par)
  over_up <- ifelse(strand(tx_par) == "+", a, b)
  over_down <- ifelse(strand(tx_par) == "+", b, a)
  weak_over_up <- weak_overlap & over_up
  weak_over_down <- weak_overlap & over_down
  mat[over_s, "Up"] <- as.logical(unlist(tapply(weak_over_up, queryHits(hits), any)))
  mat[over_s, "Down"] <- as.logical(unlist(tapply(weak_over_down, queryHits(hits), any)))
  # Switch the strand of transcripts:
  tx_sw <- tx
  strand(tx_sw) <- ifelse(strand(tx) == "+", "-", "+")
  tx_sw_start <- resize(tx_sw, 1, "end")
  # Detect DNC transcripts (start within 500 bp from TSS of a known gene):
  dnc_intervals <- suppressWarnings(trim(flank(genes, 500)))
  mat[, "Divergent"] <- tx_sw_start %over% dnc_intervals
  # Find overlaps with genes in the antisense orientation:
  over_as <- tx_sw %over% genes
  hits2 <- findOverlaps(tx_sw, genes)
  tx_sw_start_par <- tx_sw_start[queryHits(hits2)]
  genes_par2 <- genes[subjectHits(hits2)]
  # Calculate relative distance between sTSS and asTSS:
  genes_par2_start <- resize(genes_par2, 1, "start")
  dist <- width(pgap(tx_sw_start_par, genes_par2_start))
  rel_dist <- dist  / width(genes_par2)
  # Find convergent, TTS-AS and distal AS transcripts:
  conv <- rel_dist <= 0.5
  tts_as <- rel_dist > 0.5 & rel_dist <= 1.2
  distal <- rel_dist > 1.2
  mat[over_as, "Convergent"] <- as.logical(unlist(tapply(conv, queryHits(hits2), any)))
  mat[over_as, "TTS_AS"] <- as.logical(unlist(tapply(tts_as, queryHits(hits2), any)))
  mat[over_as, "Distal_AS"] <- as.logical(unlist(tapply(distal, queryHits(hits2), any)))
  # Transcripts which are FALSE in all these classes are intergenic:
  mat[, "Intergenic"] <- ifelse(rowSums(mat[, 1:7]) > 0, FALSE, TRUE)
  return(mat)
}
