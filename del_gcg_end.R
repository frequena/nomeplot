del_gcg_end <- function(refseq, discard_gcg) {
  
  raw_cg <-gregexpr("CG", refseq)
  raw_cg_n <- as.numeric(raw_cg[[1]])

  if (discard_gcg == 'nome') {
  refseq_split = strsplit(refseq , "")
  if (head(raw_cg_n, n = 1) == 1) {
    forgot_cg <- head(raw_cg_n, n = 1)
    raw_cg_n <- raw_cg_n[-1]
    for (i in raw_cg_n) {
      if (refseq_split[[1]][i-1] == "G") {
        raw_cg_n = raw_cg_n[- which(raw_cg_n == i)]
      }
    }
    raw_cg_n <- append(raw_cg_n, forgot_cg, 0)
  } else {
    for (i in raw_cg_n) {
      if (refseq_split[[1]][i-1] == "G") {
        raw_cg_n = raw_cg_n[- which(raw_cg_n == i)]
      }
    } 
  }
 }
  return(list(raw_cg_n))
}
