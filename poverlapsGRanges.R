# Unfortunately, poverlaps() function was defined in IRanges package only, but not in GenomicRanges. Thus, it cannot take GRanges objects as input;
# The function below was written to fix this issue;
# It takes two GRanges of the same length and returns a logical vector showing which ranges in "gr1" overlap with their counterparts in "gr2";

poverlapsGRanges <- function(gr1, gr2, maxgap = 0L, minoverlap = 1L, type = c("any", "start", "end", "within", "equal"), ignore.strand = FALSE) {
  require(GenomicRanges)
  stopifnot(length(gr1) == length(gr2))
  goodSeqnames <- as.character(seqnames(gr1)) == as.character(seqnames(gr2))
  if (isTRUE(ignore.strand)) {
    goodStrands <- TRUE
  } else {
    goodStrands <- GenomicRanges:::compatibleStrand(strand(gr1), strand(gr2))
  }
  goodRanges <- poverlaps(ranges(gr1), ranges(gr2), maxgap, minoverlap, type)
  return(goodSeqnames & goodStrands & goodRanges)
}
