primer_seq <- function(r_problem, fwd_primer, rv_primer, maxLmismatch = 1, maxRmismatch = 1) {
  fwd_primer <- toupper(as.character(fwd_primer))
  rv_primer <- toupper(as.character(rv_primer))
  for (i in 1:length(r_problem)) {
    result <- matchLRPatterns(Lpattern = fwd_primer, Rpattern = rv_primer, subject = r_problem[[i]],
                              max.Lmismatch = maxLmismatch, max.Rmismatch =maxRmismatch ,
                              max.gaplength = 1000)
    if (class(try(result[[1]]))  == "try-error") {
      cat(paste(names(r_problem)[i],"couldn't be trimmed","\n"))
      next
    }
    else {
      cat(names(r_problem)[i], "trimmed correctly","\n")
      r_problem[i] <- as(result[1], "DNAStringSet")
    }
  }
  return(r_problem)

}