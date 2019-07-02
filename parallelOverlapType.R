# This function takes two parallel GRanges as input and returns a factor with 7 levels:
# "up"/"down": the range in gr2 overlaps the 5'-end or 3'-end of the range in gr1;
# "inside": gr2 range is found within the gr1 range;
# "contains": gr2 range includes the gr1 range;
# "exact": gr2 and gr1 ranges are identical;
# "no_up"/"no_down": gr2 range does not overlap the gr1 range and is located upstream or downstream from the gr1 range;

parallelOverlapType <- function(gr1, gr2) {
  stopifnot(length(gr1) == length(gr2)) # gr1 and gr2 are expected to be parallel
  stopifnot(all(strand(gr1) %in% c("+", "-"))) # All intervals in gr1 are expected to have strand info. Strandness of gr2 is not taken into account.
  out <- vector("character", length(gr1))
  a <- start(gr2) <= start(gr1) & end(gr2) < end(gr1) & end(gr2) >= start(gr1) # gr2 overlaps the beginning of gr1
  b <- start(gr2) > start(gr1) & end(gr2) < end(gr1) # gr2 is within gr1
  c <- start(gr2) > start(gr1) & end(gr2) >= end(gr1) & start(gr2) <= end(gr1) # gr2 overlaps the end of gr1
  d <- start(gr2) < start(gr1) & end(gr2) > end(gr1) # gr2 includes (contains) gr1
  e <- start(gr2) == start(gr1) & end(gr2) == end(gr1)
  f <- end(gr2) < start(gr1) # no overlap
  g <- start(gr2) > end(gr1)
  no_up <- ifelse(strand(gr1) == "+", f, g)
  no_down <- ifelse(strand(gr1) == "+", g, f)
  out[no_up] <- "no_up"
  out[no_down] <- "no_down"
  up <- ifelse(strand(gr1) == "+", a, c)
  down <- ifelse(strand(gr1) == "+", c, a)
  out[up] <- "up"
  out[down] <- "down"
  out[b] <- "inside"
  out[d] <- "contains"
  out[e] <- "exact"
  return(as.factor(out))
}