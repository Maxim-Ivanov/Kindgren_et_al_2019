# This function can normalize the coverage in a GRanges object to the defined number of tags (1M by default);
# "gr" = GRanges object with mcols having the score column;
# If any intervals (typically nuclear protein-coding genes) were provided through the "by" argument, then the signal is normalized to the defined number of tags in these intervals
# (i.e. the total final number of tags may exceed the defined number);

normalizeGR <- function(gr, norm = 1000000, norm_factor = NULL, precision = 6, by = NULL) {
  gr <- gr[mcols(gr)$score > 0]
  if (is.null(norm_factor)) {
    if (!is.null(by) && class(by) == "GRanges") {
      rel_reads <- gr[gr %over% by]
    } else {
      rel_reads <- gr
    }
    auc <- sum(mcols(rel_reads)$score * width(rel_reads))
    norm_factor <- auc / norm
  } else if (!is.numeric(norm_factor) || length(norm_factor) != 1) {
    stop("Invalid norm_factor!")
  }
  mcols(gr)$score <- round(mcols(gr)$score / norm_factor, precision)
  return(gr)
}