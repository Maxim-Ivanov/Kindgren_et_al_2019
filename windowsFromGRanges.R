# Each input interval in "gr" is converted to a set of windows of fixed width with defined offset from each other;
# ID of the input interval is recorded in mcols;
# (the same effect could be achieved by tiling the whole genome and subsetting tiles which overlap with input intervals, but this is slow and memory consuming)

windowsFromGRanges <- function(gr, window_width, window_offset) {
  new_start_list <- mapply(function(x, y) {seq(from = x, to = y, by = window_offset)}, start(gr), end(gr) - window_width + 1, SIMPLIFY = FALSE)
  nums <- lengths(new_start_list)
  new_start <- unlist(new_start_list)
  new_seqnames <- Rle(values = as.factor(seqnames(gr)), lengths = nums)
  new_strand <- Rle(values = as.factor(strand(gr)), lengths = nums)
  out <- GRanges(seqnames = new_seqnames, ranges = IRanges(start = new_start, width = window_width), strand = new_strand, seqinfo = seqinfo(gr))
  mcols(out)$ID <- Rle(values = 1:length(gr), lengths = nums)
  return(out)
}
