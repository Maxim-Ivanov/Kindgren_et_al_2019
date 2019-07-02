# These functions are used to compute a metagene matrix from two GRanges objects:
# "signal" has mcols with the score column and represents NGS coverage along the chromosomes;
# "intervals" are the chosen genomic windows (e.g. 1Kb windows centered at TSS) with constant or variable width;
# Returns numeric matrix which can be used to draw either a heatmap (using pheatmap package) or a metagene plot (using drawMetagenePlot() function)


findOutOfBounds <- function(gr) {
  lookup <- as.list(seqlengths(seqinfo(gr)))
  lengths <- lookup[as.character(seqnames(gr))]
  out <- (start(gr) < 0) | (end(gr) > lengths)
  return(out)
}

calcMatrixLength <- function(intervals, matrix.length) {
  max_width <- max(width(intervals))
  min_width <- min(width(intervals))
  if (is.numeric(matrix.length)) {
    out <- matrix.length
  } else if (max_width == min_width) {
    out <- max_width
  } else if (matrix.length == "max") {
    out <- max_width
  } else if (matrix.length == "min") {
    out <- min_width
  } else if (matrix.length == "mean") {
    out <- round(mean(width(intervals)))
  } else if (matrix.length == "median") {
    out <- round(median(width(intervals)))
  } else {
    message("Not sure how to determine the number of bins. Using median interval length by default..."); flush.console()
    out <- round(median(width(intervals)))
  }
  message("Matrix length = ", out); flush.console()
  return(out)
}

expandOrShrink <- function(x, mlen) {
  if (length(x) == mlen) {
    return(as.numeric(x))
  } else {
    runLength(x) <- runLength(x) * mlen
    return(colMeans(matrix(x, ncol = mlen)))###
  }
}

trimOrFill <- function(x, mlen, anchor, na.as.zeros) {
  if (length(x) == mlen) {
    out <- x
  } else if (length(x) > mlen) {
    if (anchor == "start") {
      out <- x[1:mlen]
    } else if (anchor == "end") {
      out <- x[(length(x)-mlen+1):length(x)]
    }
  } else {
    if (isTRUE(na.as.zeros)) { vals <- 0 } else { vals <- NA }
    to_add <- Rle(values = vals, lengths = mlen - length(x))
    if (anchor == "start") {
      out <- c(x, to_add)
    } else if (anchor == "end") {
      out <- c(to_add, x)
    }
  }
  return(as.numeric(out))
}


calcMatrix <- function(signal, intervals, strand = NULL, max_w, min_w, scaling, mlen, anchor, na.as.zeros) {
  if (!is.null(strand)) {
    intervals <- intervals[strand(intervals) == strand]
    signal <- signal[strand(signal) == strand]
  }
  rlelist <- round(coverage(signal, weight = "score"), 8) # to avoid small negative and positive values instead of zeros
  rlelist_split <- rlelist[intervals]
  rlelist_split <- revElements(rlelist_split, strand(intervals) == "-")
  if ((max_w != min_w) || (max_w != mlen)) {
    if (isTRUE(scaling)) {
      if (length(rlelist_split) > 0) { message("Expanding and shrinking (this can be slow)..."); flush.console() }
      numlist <- lapply(rlelist_split, expandOrShrink, mlen = mlen)
    } else {
      if (length(rlelist_split) > 0) { message("Trimming and filling..."); flush.console() }
      numlist <- lapply(rlelist_split, trimOrFill, mlen = mlen, anchor = anchor, na.as.zeros = na.as.zeros)
    }
  } else {
    numlist <- as(rlelist_split, "NumericList")
  }
  mat <- do.call(rbind, numlist)
  return(mat)
}

average_bin <- function(x) {
  if (all(is.na(x))) {
    return(NA)
  } else {
    return(sum(x, na.rm = TRUE) / length(x))
  }
}

shrink_row <- function(x, mask) {
  return(as.numeric(tapply(x, mask, average_bin, simplify = FALSE)))
}

shrink_to_bins <- function(mat, binsize) {
  mask <- rep(seq(1, ncol(mat) / binsize), each = binsize)
  out <- t(apply(mat, 1, shrink_row, mask = mask))
  return(out)
}

metageneMatrix <- function(signal, intervals, scaling = FALSE, matrix.length = NA, anchor = "start", na.as.zeros = FALSE, 
                           skip.zeros = TRUE, skip.outliers = 0.995, skip.top.obs = FALSE, n.top.obs = 3, 
                           equal.weights = FALSE, antisenseMode = FALSE, shrink = FALSE, binsize = 5) {
  library(GenomicRanges)
  out_of_bound <- findOutOfBounds(intervals)
  if (sum(out_of_bound) > 0) {
    message(sum(out_of_bound), " intervals were out-of-bounds;"); flush.console()
    intervals <- intervals[!out_of_bound]
  }
  names(intervals) <- 1:length(intervals) # enumerate intervals
  mlen <- calcMatrixLength(intervals, matrix.length)
  max_w <- max(width(intervals)); min_w <- min(width(intervals))
  if (any(strand(intervals)=="*")) {
    message(round(sum(strand(intervals)=="*") / length(intervals) * 100, 1), "% intervals are unstranded!"); flush.console()
  }
  if (any(strand(signal)=="*")) {
    message(round(sum(strand(signal)=="*") / length(signal) * 100, 1), "% signal values are unstranded!"); flush.console()
  }
  if (isTRUE(antisenseMode)) {
    stranded <- strand(signal) %in% c("+", "-")
    strand(signal)[stranded] <- ifelse(strand(signal)[stranded]=="+", "-", "+")
    message("Antisense mode: signal was flipped to the opposite strand;"); flush.console()
  }
  interval_strands <- list("+", "-", "*")
  signal_strands <- list(c("+", "*"), c("-", "*"), c("+", "-", "*"))
  matlist <- vector("list", 3)
  for (i in 1:3) {
    curr_signal <- signal[strand(signal) %in% signal_strands[[i]]]
    curr_intervals <- intervals[strand(intervals) %in% interval_strands[[i]]]
    if (length(curr_intervals) > 0) {
      curr_mat <- calcMatrix(signal = curr_signal, intervals = curr_intervals, max_w = max_w, min_w = min_w, 
                             scaling = scaling, mlen = mlen, anchor = anchor, na.as.zeros = na.as.zeros)
      rownames(curr_mat) <- names(curr_intervals) # to trace back the original interval (because matrix rows are permuted at this step)
      matlist[[i]] <- curr_mat
    }
  }
  mat <- do.call(rbind, matlist)
  mat <- mat[order(as.numeric(rownames(mat))), ] # restore the original order of intervals
  if (isTRUE(shrink) && is.numeric(binsize) && length(binsize) == 1) {
    if (ncol(mat) %% binsize == 0) {
      message("Averaging signal within ", binsize, " bp bins;")
      mat <- shrink_to_bins(mat, binsize) # when it is required to average signal within N bp bins but scaling == FALSE
      message("Now matrix length is ", ncol(mat), "!")
    } else {
      message("Check the binsize parameter!")
    }
  }
  gene_cov <- rowSums(mat, na.rm = TRUE)
  if (is.numeric(skip.outliers) & length(skip.outliers) == 1) {
    q <- quantile(gene_cov, skip.outliers)
    outliers <- gene_cov > q
    mat <- mat[!outliers, ]
    message("Skipped ", sum(outliers), " potential outliers;"); flush.console()
    gene_cov <- rowSums(mat, na.rm = TRUE) # recalculate gene_cov
  }
  if (isTRUE(skip.top.obs) && is.numeric(n.top.obs) && length(n.top.obs) == 1) {
    mat_sorted <- apply(mat, 2, sort, decreasing = TRUE, na.last = TRUE) # sort each column independently from other columns
    max_vals <- mat_sorted[n.top.obs, ] # find N'th top value
    for (j in 1:ncol(mat)) {
      to_skip <- mat[, j] >= max_vals[j]
      mat[to_skip, j] <- NA # in each column, change N top observations to NA
    }
    message("Skipped ", n.top.obs, " top observations in each bin;")
    gene_cov <- rowSums(mat, na.rm = TRUE)
  }
  if (isTRUE(skip.zeros)) {
    zeros <- gene_cov == 0
    if (sum(!zeros) < 100) {
      stop("Too little intervals with non-zero signal! Consider changing skip.zeros to FALSE!")
    }
    if (sum(zeros) > 0) {
      mat <- mat[!zeros, ]
      message("Skipped ", sum(zeros), " intervals with zero signal;"); flush.console()
    }
  }
  if (isTRUE(equal.weights)) {
    mat <- t(apply(mat, 1, function(x) { x / sum(x, na.rm=TRUE) }))
    message("All intervals were assigned equal weights;"); flush.console()
  }
  return(mat)
}
