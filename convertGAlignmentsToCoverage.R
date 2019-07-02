# This script basically does the same job as BEDtools genomecov;
# "ga" = GAlignments object produced from a single-end BAM file by GenomicAlignments::readGAlignments()
# Returns GRanges object;

convertGAlignmentsToCoverage <- function(ga, split = TRUE, mode = "whole_read", merge.strands = FALSE, flip.strands = FALSE) {
  require(GenomicAlignments)
  if (mode %in% c("start", "end")) {
    if (class(ga)=="GAlignments") {
      gr <- granges(ga)
    } else if (class(ga)=="GAlignmentPairs") {
      gr <- granges(c(first(ga), last(ga)))
    }
    gr <- resize(gr, width=1, fix=mode)
  } else if (mode=="whole_read") {
    if (isTRUE(split)) {
      gr <- unlist(grglist(ga))
    } else {
      gr <- granges(ga)
    }
  }
  if (isTRUE(merge.strands)) {
    cov_all <- coverage(gr)
    out <- bindAsGRanges(score=cov_all)
  } else {
    cov_fw <- coverage(gr[strand(gr)=="+"])
    cov_rev <- coverage(gr[strand(gr)=="-"])
    out_fw <- bindAsGRanges(score=cov_fw)
    strand(out_fw) <- "+"
    out_rev <- bindAsGRanges(score=cov_rev)
    strand(out_rev) <- "-"
    out <- c(out_fw, out_rev)
    if (isTRUE(flip.strands)) {
      strand(out) <- ifelse(strand(out)=="+", "-", "+")
    }
  }
  out <- out[score(out)>0]
  out <- sortSeqlevels(out)
  out <- out[order(seqnames(out), start(out))]
  return(out)
}