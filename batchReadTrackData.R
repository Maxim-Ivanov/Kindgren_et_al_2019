# This function iterates over a list of Bedgraph or Bigwig filenames and calls rtracklayer::import() function for each filename;
# A GRangesList object is returned;
# This function can handle different ways of encoding strandness:
# a) Bedgraph files may contain strand info encoded as positive and negative scores in the 4th column;
# b) If the files are either Bedgraphs with only positive scores or BigWigs (which cannot contain negative scores), then they are expected to have the same strandness (unknown by default);
# In this case, the strand info has to be explicitly provided through the "strand" argument;
# For Bedgraph files, it is always highly recommended to provide a valid Seqinfo object through the "seqinfo" argument.

batchReadTrackData <- function(filenames, dir = ".", format = NULL, strand = "*", seqinfo = NULL) {
  require(rtracklayer)
  output <- vector("list", length(filenames))
  for (i in seq_along(filenames)) {
    fn <- filenames[[i]]
    cat(fn, "\n"); flush.console()
    if (is.null(format)) {
      if (is.null(seqinfo)) {
        gr <- trim(import(file.path(dir, fn)))
      } else {
        gr <- trim(import(file.path(dir, fn), seqinfo=seqinfo))
      }
    } else {
      if (is.null(seqinfo)) {
        gr <- trim(import(file.path(dir, fn), format=format))
      } else {
        gr <- trim(import(file.path(dir, fn), format=format, seqinfo=seqinfo))
      }
    }
    if (any(score(gr) < 0)) {
      strand(gr) <- ifelse(score(gr) >= 0, "+", "-")
      score(gr) <- abs(score(gr))
    } else {
      strand(gr) <- strand
    }
    output[[i]] <- sort(gr)
    names(output)[[i]] <- basename(fn)
  }
  output <- GRangesList(output)
  return(output)
}
