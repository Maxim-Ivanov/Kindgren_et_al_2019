saveGRangesAsBedGraph <- function(gr, filepath, colname = "score") {
  library(rtracklayer)
  mcols(gr)[[colname]] <- ifelse(strand(gr)=="-", -mcols(gr)[[colname]], mcols(gr)[[colname]])
  strand(gr) <- "*"
  gr <- sort(gr)
  con <- gzfile(filepath, "w")
  writeLines("track type=bedGraph color=0,0,0 altColor=128,128,128", con)
  rtracklayer::export.bedGraph(gr, con)
  close(con)
}