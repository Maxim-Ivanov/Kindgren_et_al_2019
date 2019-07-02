# This function is used to post-process plaNET-Seq and pNET-Seq BAM files returned by 01-Alignment_plaNET-Seq.sh

loadNETSeqBAM <- function(bamfiles, bamdir, ss_5p = NULL, ss_3p = NULL, bs = NULL, pattern = "TGGAATTCTC", trim_names = NULL, expand_names = NULL, 
                             remove_split_reads = TRUE, remove_ss_reads = TRUE, remove_misprimed = TRUE, save_rds = FALSE) {
  library(GenomicAlignments)
  mat <- matrix(0, nrow = length(bamfiles), ncol = 6)
  for (i in seq_along(bamfiles)) {
    res <- vector("integer", 6)
    bamfile <- bamfiles[[i]]
    message(bamfile); flush.console()
    data <- readGAlignments(file.path(bamdir, bamfile), param = ScanBamParam(mapqFilter = 10))
    message("\t", length(data), " reads loaded;")
    res[[1]] <- length(data)
    strand(data) <- ifelse(strand(data)=="+", "-", "+") # flip the strand info (because this is NET-Seq)
    if (isTRUE(remove_split_reads)) {
      split_reads <- njunc(data) > 0
      message("\t", sum(split_reads), " of them are split;")
      res[[2]] <- sum(split_reads)
      data <- data[!split_reads]
    } else {
      res[[2]] <- 0
    }
    # Filter out reads which end on splice sites:
    if (isTRUE(remove_ss_reads) & !is.null(ss_5p) & !is.null(ss_3p)) {
      pol_pos <- resize(granges(data), 1, "end")
      over_ss_3p <- overlapsAny(pol_pos, ss_3p)
      message("\t", sum(over_ss_3p), " reads end in 3p splice sites;")
      res[[3]] <- sum(over_ss_3p)
      over_ss_5p <- overlapsAny(pol_pos, ss_5p)
      message("\t", sum(over_ss_5p), " reads end in 5p splice sites;")
      res[[4]] <- sum(over_ss_5p)
      data <- data[!over_ss_5p & !over_ss_3p]
    } else {
      res[[3]] <- 0
      res[[4]] <- 0
    }
    
    # Filter out misprimed reads:
    if (isTRUE(remove_misprimed) & !is.null(bs) && class(bs) == "BSgenome") {
      # Get 8 bp immediately downstream the end of flipped reads:
      downstr <- suppressWarnings(trim(flank(granges(data), 10, start = FALSE)))
      seq <- as(Views(bs, downstr), "DNAStringSet")
      misprimed <- vcountPattern(pattern, seq, max.mismatch = 2) == 1
      message("\t", sum(misprimed), " reads were suspected to be misprimed;")
      res[[5]] <- sum(misprimed)
      data <- data[!misprimed]
    } else {
      res[[5]] <- 0
    }
    # Produce output filename:
    out_filename <- sub(".bam$", "", bamfile)
    if (!is.null(trim_names) & is.character(trim_names)) {
      for (j in seq_along(trim_names)) {
        out_filename <- sub(trim_names[[j]], "", out_filename)
      }
    }
    if (!is.null(expand_names) & is.character(expand_names)) {
      out_filename <- paste0(out_filename, expand_names)
    }
    # Save GAlignments as RDS:
    if (isTRUE(save_rds)) {
      saveRDS(data, paste0(out_filename, ".RDS"))
    }
    # Save NET-Seq coverage as Bedgraph:
    cov <- convertGAlignmentsToCoverage(data, mode = "end")
    saveGRangesAsBedGraph(cov, paste0(out_filename, ".bedgraph.gz"))
    message("\t", length(data), " clean reads returned;")
    res[[6]] <- length(data)
    mat[i, ] <- res
  }
  df <- as.data.frame(mat)
  names(df) <- c("Loaded", "Split", "3pSS", "5pSS", "Misprimed", "Returned")
  row.names(df) <- sub(".bam", "", bamfiles)
  return(df)
}
