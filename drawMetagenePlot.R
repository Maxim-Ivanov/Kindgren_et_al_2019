# These scripts are required to draw a metagene plot from a list of metagene matrices produced by metageneMatrix()

smoothVector <- function(vec, method, degree) {
  if (method == "tukey") {
    out <- smooth(vec)
  } else if (method == "lowess") {
    out <- lowess(1:length(vec), vec, f = degree)$y
  } else if (method == "friedman") {
    out <- supsmu(1:length(vec), vec, bass = degree)
  }
  return(out)
}

drawMetagenePlot <- function(mat_list, drawCI=TRUE, x.axis=FALSE, title = "Metagene_plot", filename = NULL, xlabel = "Intervals of interest", 
                             ylabel = "Average signal", vline = FALSE, hline = FALSE, width=NA, height=NA, units="cm", plotPDF=TRUE, scale.plot=NA, linetype=0, 
                             alpha=0.25, ylim=NA, smoothing = FALSE, method = "tukey", degree = NA, custom.colors = NA, out_dir = ".") {
  library(ggplot2)
  if (identical(x.axis, FALSE)) {
    x.axis <- seq_len(ncol(mat_list[[1]]))
  }
  stopifnot(length(x.axis) == ncol(mat_list[[1]]))
  long_df <- data.frame()
  for (i in seq(length(mat_list))) {
    curr_mat <- mat_list[[i]]
    curr_name <- names(mat_list)[[i]]
    if (!is.na(scale.plot) && is.numeric(scale.plot) && length(scale.plot)==length(mat_list)) {
      curr_mat <- curr_mat * scale.plot[[i]]
    }
    avg <- apply(curr_mat, 2, mean, na.rm = TRUE)
    message("Max value ", i, " = ", max(avg, na.rm = TRUE)); flush.console()
    if (isTRUE(smoothing)) {
      avg <- smoothVector(avg, method = method, degree = degree)
    }
    df <- data.frame("pos" = x.axis, "avg" = avg, "group" = curr_name)
    if (isTRUE(drawCI)) {
      sem <- apply(curr_mat, 2, sd, na.rm = TRUE) / sqrt(nrow(curr_mat))
      if (isTRUE(smoothing)) {
        sem <- smoothVector(sem, method = method, degree = degree)
      }
      lower <- avg - 1.96 * sem; upper <- avg + 1.96 * sem
      df <- cbind(df, data.frame("lower" = lower, "upper" = upper))
    }
    long_df <- rbind(long_df, df)
  }
  p <- ggplot(long_df, aes(x = pos, y = avg, color = group)) + geom_line(size=1) + ggtitle(title) + xlab(xlabel) + ylab(ylabel) + theme_bw()
  if (isTRUE(drawCI)) {
    p <- p + geom_ribbon(aes(x = pos, ymin = lower, ymax = upper, fill = group), alpha = alpha, linetype = linetype)
  }
  if (is.numeric(vline)) {
    p <- p + geom_vline(xintercept = vline)
  }
  if (is.numeric(hline)) {
    p <- p + geom_hline(yintercept = hline)
  }
  if (!all(is.na(ylim))) {
    p <- p + ylim(ylim[[1]], ylim[[2]])
  }
  if (!all(is.na(custom.colors))) {
    p <- p + scale_colour_manual(values = custom.colors) + scale_fill_manual(values = custom.colors)
  }
  if (is.null(filename)) {
    filename <- title
  }
  filename <- sub("\n", " ", filename)
  ggsave(filename = file.path(out_dir, paste0(filename, ".png")), plot = p, width=width, height=height, units=units)
  if (isTRUE(plotPDF)) {
    ggsave(filename = file.path(out_dir, paste0(filename, ".pdf")), plot = p, width=width, height=height, units=units)
  }
}

