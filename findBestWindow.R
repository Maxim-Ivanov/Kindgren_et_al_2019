# This function is used when calculating TSS stalling index (TSS SI) and intronic stalling index (ISI);

findBestWindow <- function(df, random = FALSE) {
  max_cov <- df$cov == max(df$cov) # matrix instead of vector...
  if (sum(max_cov) == 1) {
    best_idx <- df$idx[max_cov]
  } else {
    if (isTRUE(random)) {
      best_idx <- df$idx[sample(which(max_cov), size = 1)]
    } else {
      dists <- df$dist[max_cov]
      min_dist <- abs(dists) == min(abs(dists)) # fixed bug on 2019-06-17
      if (sum(min_dist) > 1) { # fixed bug on 2019-06-17
        min_dist <- sample(which(min_dist), size = 1)
      }
      best_idx <- df$idx[max_cov][min_dist]
    }
  }
  return(best_idx)
}