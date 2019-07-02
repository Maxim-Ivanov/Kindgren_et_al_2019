# These functions are required to call the read-through distance by 08-Readthrough_distance.R

find_end_of_TRUE_run <- function(x, mode) {
  if (mode == "left_to_right") {
    if (all(x)) {
      out <- length(x)
    } else if (all(!x)) {
      out <- NA
    } else {
      out <- min(which(!x)) - 1
      if (out == 0) { out <- NA }
    }
  } else if (mode == "right_to_left") {
    if (all(x)) {
      out <- NA
    } else if (all(!x)) {
      out <- length(x)
    } else {
      out <- max(which(!x)) + 1
      if (out > length(x)) { out <- NA }
    }
  }
  return(out)
}

smoothed_ecdf <- function(vals, adjust = 1) {
  dens <- density(vals, adjust = adjust)
  dens$y <- cumsum(dens$y) / sum(dens$y) # convert PDF to CDF
  return(dens)
}

computeRTtails <- function(genes, cov_mcol, win_cov, lambda_ctrl) {
  # Extract windows covering the genes:
  hits <- findOverlaps(win_cov, genes, type = "within") # the reverse order of query and subject because of type = "within"
  hits <- Hits(from = subjectHits(hits), to = queryHits(hits), nLnode = nRnode(hits), nRnode = nLnode(hits)) # flip queryHits() and subjectHits() for convenience
  hits <- hits[order(queryHits(hits))] # sort by queryHits
  # Mask genes which do not contain sufficient number of genic windows:
  rle <- as(queryHits(hits), "Rle")
  runlen <- runLength(rle)
  runval <- runValue(rle)
  bad_gene_idx <- runval[runlen < 10]
  zero_gene_idx <- setdiff(1:length(genes), unique(queryHits(hits))) # genes with no windows at all
  mask <- rep(TRUE, length(genes))
  mask[c(bad_gene_idx, zero_gene_idx)] <- FALSE # TRUE for good genes
  hits <- hits[!queryHits(hits) %in% bad_gene_idx]
  # Compute the distribution of tag counts for each gene:
  cov_genic <- mcols(win_cov)[subjectHits(hits), cov_mcol]
  cdf <- tapply(cov_genic, queryHits(hits), smoothed_ecdf, simplify = FALSE)
  # Extract windows within the gaps between gene ends/PAS and downstream TSS:
  hits2 <- findOverlaps(win_cov, mcols(genes[mask])$gap, type = "within")
  hits2 <- Hits(from = subjectHits(hits2), to = queryHits(hits2), nLnode = nRnode(hits2), nRnode = nLnode(hits2))
  hits2 <- hits2[order(queryHits(hits2))] # (queryHits point to gene indexes in the masked subset)
  # For each window, compute the likelihood under the gene model:
  cov_gap <- mcols(win_cov)[subjectHits(hits2), cov_mcol]
  cov_gap_spl <- split(cov_gap, queryHits(hits2))
  prob_genic <- unname(unlist(mapply(function(a, b) { return(approx(a$x, a$y, xout = b, yleft = 0, yright = 1)$y) }, cdf, cov_gap_spl, SIMPLIFY = FALSE)))
  # Compute also the likelihood under the control (random) model:
  prob_ctrl <- ppois(cov_gap, lambda = lambda_ctrl, lower.tail = FALSE)
  # Calculate the ratio of probabilities:
  ratio <- prob_genic / prob_ctrl
  ratio[is.infinite(ratio)] <- 1
  ratio[is.nan(ratio)] <- 1
  # Make decisions:
  decisions <- relist(ratio >= 1, cov_gap_spl)
  # Find position of the TRUE immediately before the leftmost FALSE:
  idx1 <- unname(unlist(lapply(decisions, find_end_of_TRUE_run, mode = "left_to_right")))
  # Alternatively, find the TRUE immediately after the rightmost FALSE:
  idx2 <- unname(unlist(lapply(decisions, find_end_of_TRUE_run, mode = "right_to_left")))
  # Choose results based on strand orientation of the original genes:
  idx <- ifelse(strand(genes[mask]) == "+", idx1, idx2)
  # Additionally mask genes with idx == NA (no RT tail at all):
  not_na <- !is.na(idx)
  idx <- idx[not_na]
  mask[mask][!not_na] <- FALSE
  hits2 <- hits2[queryHits(hits2) %in% which(not_na)]
  # Translate the local indexes to the global ones:
  rle2 <- as(queryHits(hits2), "Rle")
  runlen2 <- runLength(rle2)
  base <- cumsum(c(0, runlen2[-length(runlen2)]))
  global_idx <- base + idx
  # Extract the last windows in RT tails:
  win <- win_cov[subjectHits(hits2)]
  win_last <- win[global_idx]
  # Insert last windows into mcols(genes):
  dummy_tail <- resize(mcols(genes)$gap, 1, "start")
  real_tail <- resize(win_last, 1, "end") # length(real_tail) < length(dummy_tail)
  coords <- start(dummy_tail)
  coords[mask] <- start(real_tail)
  rt_tail_ends <- GRanges(seqnames = seqnames(genes), ranges = IRanges(start = coords, width = 1), strand = strand(genes))
  rt_tails <- pgap(resize(mcols(genes)$gap, 1, "start"), rt_tail_ends)
  mcols(genes)[paste0("RT_", cov_mcol)] <- rt_tails
  return(genes)
}