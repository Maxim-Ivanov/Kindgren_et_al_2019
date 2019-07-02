# "ebg" = GRangesList containing coordinates of exons grouped by gene: ebg <- exonsBy(txdb, by = "gene")
# Returns list of length 3:
# The first element is GRanges object containing all exons except for the first ones;
# The second element is GRanges with all exons except for the last ones;
# The third element is GRanges with all exons except for both first and last ones;
# This script considers continuous exonic intervals belonging to the same gene. Thus, a gene may contain multiple first and/or last exons. All of them are detected and skipped;

removeFirstAndLastExons <- function(ebg) {
  ebg <- sort(ebg)
  ebg_unl <- unlist(ebg)
  ebg_red <- reduce(ebg)
  ebg_red_unl <- unlist(ebg_red)
  enum_fw <- unlist(sapply(lengths(ebg_red), function(x) { seq(1, x) }))
  enum_rev <- unlist(sapply(lengths(ebg_red), function(x) { seq(x, 1) }))
  first_red <- ebg_red_unl[as.logical((strand(ebg_red_unl) == "+" & enum_fw == 1) | (strand(ebg_red_unl) == "-" & enum_rev == 1))]
  h1 <- findOverlaps(ebg_unl, first_red)
  same_gene_1 <- names(ebg_unl[queryHits(h1)]) == names(first_red[subjectHits(h1)])
  idx_1 <- unique(queryHits(h1[same_gene_1]))
  exons_wo_first <- ebg_unl[-idx_1]
  last_red <- ebg_red_unl[as.logical((strand(ebg_red_unl) == "+" & enum_rev == 1) | (strand(ebg_red_unl) == "-" & enum_fw == 1))]
  h2 <- findOverlaps(ebg_unl, last_red)
  same_gene_2 <- names(ebg_unl[queryHits(h2)]) == names(last_red[subjectHits(h2)])
  idx_2 <- unique(queryHits(h2[same_gene_2]))
  exons_wo_last <- ebg_unl[-idx_2]
  internal_exons <- ebg_unl[-unique(c(idx_1, idx_2))]
  return(list(exons_wo_first, exons_wo_last, internal_exons))
}

