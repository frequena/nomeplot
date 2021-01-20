
del_gcg_exo <- function(refseq) {
  
  raw_gc <-gregexpr("GC", refseq)
  raw_gc[[1]] <- raw_gc[[1]] + 1
  raw_gc_n <- as.numeric(raw_gc[[1]])
  last_gc <- tail(raw_gc_n, n = 1)
  refseq_split <- strsplit(paste(refseq) , "")
  
  if (tail(raw_gc_n, n = 1) == length(refseq_split[[1]])) {
    forgot_gc <- raw_gc_n[length(raw_gc_n)]
    raw_gc_n <- raw_gc_n[-length(raw_gc_n)]
    for (i in raw_gc_n) {
      if (refseq_split[[1]][i+1] == "G") {
        raw_gc_n <- raw_gc_n[- which(raw_gc_n == i)]
      }
    }
    raw_gc_n <- append(raw_gc_n,forgot_gc )
  } else {
    
  for (i in raw_gc_n) {
    if (refseq_split[[1]][i+1] == "G") {
      raw_gc_n <- raw_gc_n[- which(raw_gc_n == i)]
    }
   }
 }
  return(list(raw_gc_n))
}