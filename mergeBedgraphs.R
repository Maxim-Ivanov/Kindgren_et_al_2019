# These functions allow to combine the coverage in multiple Bedgraph files into a single Bedgraph;
# This is useful when merging technical or biological replicates;

findOutOfBounds <- function(gr) {
  lookup <- as.list(seqlengths(seqinfo(gr)))
  lengths <- lookup[as.character(seqnames(gr))]
  out <- (start(gr) < 0) | (end(gr) > lengths)
  return(out)
}

loadBgCov <- function(bg, seqinfo, precision) {
  library(rtracklayer)
  gr <- rtracklayer::import.bedGraph(bg, seqinfo=seqinfo)
  if (any(mcols(gr)$score < 0)) {
    gr_pos <- gr[mcols(gr)$score >= 0]
    gr_neg <- gr[mcols(gr)$score < 0]
    mcols(gr_neg)$score <- abs(mcols(gr_neg)$score)
    return(list(round(coverage(gr_pos, weight = "score"), precision), round(coverage(gr_neg, weight = "score"), precision)))
  } else {
    return(list(round(coverage(gr, weight = "score"), precision)))
  }
}

mergeBedgraphs <- function(bg_fw = list(), bg_rev = list(), bg_fw_rev = list(), mode = "stranded", seqinfo, norm = 1000000, precision = 8) {
  if (mode == "stranded") {
    cov_fw <- vector("list", length(bg_fw) + length(bg_fw_rev))
    cov_rev <- vector("list", length(bg_rev) + length(bg_fw_rev))
    i <- j <- 0
    if (length(bg_fw) > 0) {
      for (i in seq_along(bg_fw)) {
        cov_list <- loadBgCov(bg_fw[[i]], seqinfo = seqinfo, precision = precision)
        cov_fw[[i]] <- cov_list[[1]]
      }
    }
    if (length(bg_rev) > 0) {
      for (j in seq_along(bg_rev)) {
        cov_list <- loadBgCov(bg_rev[[j]], seqinfo = seqinfo, precision = precision)
        cov_rev[[j]] <- cov_list[[1]]
      }
    }
    if (length(bg_fw_rev) > 0) {
      for (k in seq_along(bg_fw_rev)) {
        cov_list <- loadBgCov(bg_fw_rev[[k]], seqinfo = seqinfo, precision = precision)
        cov_fw[[i + k]] <- cov_list[[1]]
        cov_rev[[j + k]] <- cov_list[[2]]
      }
    }
    cov_fw_merged <- Reduce(`+`, cov_fw)
    cov_rev_merged <- Reduce(`+`, cov_rev)
    gr_fw <- bindAsGRanges(score = cov_fw_merged)
    strand(gr_fw) <- "+"
    gr_rev <- bindAsGRanges(score = cov_rev_merged)
    strand(gr_rev) <- "-"
    gr <- c(gr_fw, gr_rev)
  } else if (mode == "unstranded") {
    all_bg <- c(bg_fw, bg_rev, bg_fw_rev)
    cov <- vector("list", length(all_bg))
    if (length(all_bg) > 0) {
      for (i in seq_along(all_bg)) {
        cov_list <- loadBgCov(bg_fw[[i]], seqinfo = seqinfo, precision = precision)
        if (length(cov_list) == 1) {
          cov[[i]] <- cov_list[[1]]
        } else {
          cov[[i]] <- c(cov_list[[1]], cov_list[[2]])
        }
      }
    }
    cov_merged <- Reduce(`+`, cov)
    gr <- bindAsGRanges(score = cov_merged)
    strand(gr) <- "*"
  }
  out <- findOutOfBounds(gr)
  if (any(out)) {
    gr <- gr[!out]
  }
  gr <- gr[order(seqnames(gr), start(gr))]
  if (is.numeric(norm) && length(norm)==1) {
    total <- round(sum(mcols(gr)$score * width(gr)))
    message(total, " tags in total;"); flush.console()
    norm_f <- norm / total
    message("Normalization factor = ", round(norm_f, precision), ";"); flush.console()
    mcols(gr)$score <- round(mcols(gr)$score * norm_f, precision)
  }
  gr <- gr[mcols(gr)$score > 0] # suppress intervals with zero scores
  return(gr)
}
