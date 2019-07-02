# This function takes GRanges object as input and returns also a GRanges object;
# The input GRanges is allowed to have ranges with width above 1, the output GRanges may only contain ranges with unit width;
# Each input ranges with width above 1 gets split into multiple ranges having width 1 and the same score;
# The unit width GRanges can be further converted to BigWig and used for as input for the CAGEfightR package;
# (due to unclear reason, CAGEfightR cannot process input files containing ranges with width above 1);

expandGRtoUnitWidth <- function(gr) {
  exp_seqnames <- rep(seqnames(gr), width(gr))
  exp_start <- unlist(mapply(function(x, y) { seq(x, x + y - 1) }, start(gr), width(gr), SIMPLIFY = FALSE))
  exp_strand <- rep(strand(gr), width(gr))
  exp_score <- rep(score(gr), width(gr))
  out <- GRanges(seqnames = exp_seqnames, ranges = IRanges(start = exp_start, width = 1), 
                 strand = exp_strand, score = exp_score, seqinfo = seqinfo(gr))
  return(out)
}
