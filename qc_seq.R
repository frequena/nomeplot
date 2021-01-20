qc_seq <- function(r_master, r_problem, convert = NULL) {
  require(Biostrings)
  data_qc <- data.frame("original" = integer(),
                        "complement" = integer(), "reverse" = integer(), "comp_reverse" = integer(),
                        "Orientation" = character(), stringsAsFactors = FALSE)
 
  mat <- matrix(0, length(DNA_ALPHABET), length(DNA_ALPHABET))
  mat[1:4, 1:4] <- c(1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1)
  rownames(mat) <- colnames(mat) <- DNA_ALPHABET[1:length(DNA_ALPHABET)]
  
  for (i in 1:length(r_problem)) {
    Temp_problem <- r_problem[[i]]
    o <- pairwiseAlignment(Temp_problem ,r_master , type = "local-global"
                           , substitutionMatrix = mat)
    c <- pairwiseAlignment(complement(Temp_problem) ,r_master , type = "local-global"
                           , substitutionMatrix = mat)
    r <- pairwiseAlignment(reverse(Temp_problem) ,r_master , type = "local-global"
                          , substitutionMatrix = mat)
    c_r <-pairwiseAlignment(reverseComplement(Temp_problem) ,r_master , type = "local-global"
                             , substitutionMatrix = mat)
    qc <- c(o@score, c@score, r@score, c_r@score)
    data_qc[i,][,1:4] <- qc
    if (max(as.numeric(data_qc[i,][,1:4])) == as.numeric(data_qc[i,][,1])) {
      data_qc[i,][,5] <- "ALigned correctly"
    }
    else {
      if (max(as.numeric(data_qc[i,][,1:4])) == as.numeric(data_qc[i,][,2])) {
          data_qc[i,][,5] <- "Maybe it's a complement sequence. Check out!"    
          if (!is.null(convert)) {
            r_problem[[i]] <- complement(Temp_problem)
          }
      }
      if (max(as.numeric(data_qc[i,][,1:4])) == as.numeric(data_qc[i,][,3])) {
          data_qc[i,][,5] <- "Maybe it's a reverse sequence. Check out!"  
          if (!is.null(convert)) {
            r_problem[[i]] <- reverse(Temp_problem)
          }
      }
      if (max(as.numeric(data_qc[i,][,1:4])) == as.numeric(data_qc[i,][,4])) {
        data_qc[i,][,5] <- "Maybe it's a complement reverse sequence. Check out!"
        if (!is.null(convert)) {
          r_problem[[i]] <- DNAString(reverseComplement(Temp_problem))
        }
      }
    }
  }
  rownames(data_qc) <- r_problem@ranges@NAMES
  result <- list("QC" = data_qc, "Problem sequences" = r_problem)
  return(result)
}


