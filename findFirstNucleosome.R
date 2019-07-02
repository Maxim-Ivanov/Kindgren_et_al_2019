# "nucl" = GRanges with all nucleosome positions;
# "genes" = GRanges containing the gene model;
# Returns a list of two GRanges objects:
# The first GRanges contains coordinates of the first nucleosomes;
# The second GRanges contains coordinates of the relevant genes (it may be shorter than the input "genes")

findFirstNucleosome <- function(nucl, genes, max_dist = 500) {
  # Center of the nucleosome has to be found within 500 bp downstream from the TSS;
  tss <- resize(granges(genes), 1, "start")
  nucl_c <- resize(nucl, 1, "center")
  idx <- precede(tss, nucl_c)
  tss <- tss[!is.na(idx)]
  genes <- genes[!is.na(idx)]
  idx <- idx[!is.na(idx)]
  nucl_c_par <- nucl[idx]
  strand(nucl_c_par) <- strand(tss)
  good <- width(pgap(tss, nucl_c_par)) <= max_dist
  genes <- genes[good]
  idx <- idx[good]
  first <- nucl[idx]
  strand(first) <- strand(genes) # assign strandness
  return(list(first, genes))
}