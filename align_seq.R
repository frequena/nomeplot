
align_seq <- function(r_master, r_problem, meth_end = NULL, meth_exo = NULL, discard_gcg) {
  
  mat <- matrix(0, length(DNA_ALPHABET), length(DNA_ALPHABET))
  mat[1:4, 1:4] <- c(1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 
                     0, 0, 1)
  rownames(mat) <- colnames(mat) <- DNA_ALPHABET[1:length(DNA_ALPHABET)]
  
  mult_align <- DNAMultipleAlignment()
  
  r_master <<- r_master
  r_problem <<- r_problem
  
  for (i in 1:length(r_problem)) {
    temp_align <- pairwiseAlignment(r_problem[i] , r_master, substitutionMatrix = mat, 
                                    type = "local-global")
    mult_align <- append(mult_align, temp_align)
  }
  mult_align <- rev(mult_align)
  mult_align[length(r_problem) + 1] <- NULL

  if (!is.null(meth_end)) {
    
    positionCG <- del_gcg_end(paste(r_master), discard_gcg)
    numberCG <- length(positionCG[[1]])
    methCG <- matrix(nrow = length(r_problem), ncol = numberCG)
    
    for (i in 1:length(r_problem)) {
      methCGTemp <- cgMethFinder_mod(r_master, Biostrings::pattern(mult_align[[i]]), discard_gcg)
      methCG[i,] <- methCGTemp
    }    
    df_end <- data.frame(methCG, stringsAsFactors = FALSE)
    colnames(df_end) <- positionCG[[1]]
    #df_end <- clean_unaligned(mult_align, df_end)
  }
  if (!is.null(meth_exo)) {
    positionGC <- del_gcg_exo(paste(r_master))
    numberGC <- length(positionGC[[1]])
    methGC <- matrix(nrow = length(r_problem), ncol = numberGC)
    
    for (i in 1:length(r_problem)) {
      methGCTemp <- cgMethFinder_exo(r_master, Biostrings::pattern(mult_align[[i]]))
      methGC[i,] <- methGCTemp
    }    
    df_exo <- data.frame(methGC, stringsAsFactors = FALSE)
    colnames(df_exo) <- positionGC[[1]]
  }
  result <- list("MethCG" = df_end, "MethGC" = df_exo)
  return(result)
}
