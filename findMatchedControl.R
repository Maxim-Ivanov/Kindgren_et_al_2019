# Both "treat" and "ctrll" are vectors of numeric values (e.g. gene expression values);
# The script tries to pick up a subset of "ctrl" values which are distributed in (more or less) the same way as the the whole set of "treat" values;
# Returns integer vector of indexes for the chosen "ctrl" values;

findMatchedControl <- function(treat, ctrl, log = FALSE, pseudocount = 0.1, nbins = 100, multiply = 1) {
  if (isTRUE(log)) {
    treat <- log2(treat + pseudocount)
    ctrl <- log2(ctrl + pseudocount)
  }
  pooled <- c(treat, ctrl)
  breaks <- seq(min(pooled), max(pooled), length.out = nbins)
  treat_cut <- cut(treat, breaks = breaks, include.lowest = TRUE)
  ctrl_cut <- cut(ctrl, breaks = breaks, include.lowest = TRUE)
  df_treat <- as.data.frame(table(treat_cut))
  names(df_treat) <- c("bin", "freq")
  df_treat$bin <- as.numeric(df_treat$bin)
  df_treat <- df_treat[df_treat$freq > 0, ]
  df_treat$freq <- df_treat$freq * multiply
  df_ctrl <- data.frame("value" = ctrl, "index" = 1:length(ctrl), "bin" = as.numeric(ctrl_cut))
  df_ctrl$bin <- as.numeric(df_ctrl$bin)
  df_ctrl <- df_ctrl[!is.na(df_ctrl$bin), ]
  results <- vector("list", nrow(df_treat))
  for (i in 1:nrow(df_treat)) {
    freq <- df_treat$freq[i]
    bin <- df_treat$bin[i]
    indexes <- df_ctrl$index[df_ctrl$bin == bin]
    if (length(indexes) == 0) {
      results[[i]] <- rep(NA, length.out = freq)
    } else if (length(indexes) == 1) {
      results[[i]] <- rep(indexes, length.out = freq)
    } else if (freq <= length(indexes)) {
      results[[i]] <- sample(indexes, size = freq)
    } else {
      results[[i]] <- sample(indexes, size = freq, replace = TRUE)
    }
  }
  out <- sort(do.call(c, results))
  return(out)
}