
collect_seq <- function(problem_seq, path_file, start_interval, end_interval, chr_loc) {
  
  
  
  flag147 <- scanBamFlag(isPaired = TRUE, isProperPair = TRUE, isUnmappedQuery = FALSE,
                         hasUnmappedMate = FALSE, isMinusStrand = TRUE, isMateMinusStrand = FALSE,
                         isFirstMateRead = FALSE, isSecondMateRead = TRUE, isNotPrimaryRead = FALSE,
                         isSecondaryAlignment = FALSE, isNotPassingQualityControls = FALSE,
                         isDuplicate = FALSE)
  
  flag99 <- scanBamFlag(isPaired = TRUE, isProperPair = TRUE, isUnmappedQuery = FALSE,
                        hasUnmappedMate = FALSE, isMinusStrand = FALSE, isMateMinusStrand = TRUE,
                        isFirstMateRead = TRUE, isSecondMateRead = FALSE, isNotPrimaryRead = FALSE,
                        isSecondaryAlignment = FALSE, isNotPassingQualityControls = FALSE,
                        isDuplicate = FALSE)
  
  flag83 <- scanBamFlag(isPaired = TRUE, isProperPair = TRUE, isUnmappedQuery = FALSE,
                        hasUnmappedMate = FALSE, isMinusStrand = TRUE, isMateMinusStrand = FALSE,
                        isFirstMateRead = TRUE, isSecondMateRead = FALSE, isNotPrimaryRead = FALSE,
                        isSecondaryAlignment = FALSE, isNotPassingQualityControls = FALSE,
                        isDuplicate = FALSE)
  
  flag163 <- scanBamFlag(isPaired = TRUE, isProperPair = TRUE, isUnmappedQuery = FALSE,
                         hasUnmappedMate = FALSE, isMinusStrand = FALSE, isMateMinusStrand = TRUE,
                         isFirstMateRead = FALSE, isSecondMateRead = TRUE, isNotPrimaryRead = FALSE,
                         isSecondaryAlignment = FALSE, isNotPassingQualityControls = FALSE,
                         isDuplicate = FALSE)
  
  flag16 <- scanBamFlag(isPaired = FALSE, isProperPair = FALSE, isUnmappedQuery = FALSE,
                         hasUnmappedMate = FALSE, isMinusStrand = TRUE, isMateMinusStrand = FALSE,
                         isFirstMateRead = FALSE, isSecondMateRead = FALSE, isNotPrimaryRead = FALSE,
                         isSecondaryAlignment = FALSE, isNotPassingQualityControls = FALSE,
                         isDuplicate = FALSE)
  
  
  
  
  # CUIDADO CUANDO HAY 0 DE ALGÚN TIPO DE SECUENCIA
  bamPath <- path_file
  
  
  what <- c("rname", "strand", "pos", "qwidth", "seq")
  which <- GRanges( seqnames = paste0('chr', chr_loc), ranges =  IRanges(as.numeric(start_interval), as.numeric(end_interval)))
  
  # FLAG 147
  
  param <- ScanBamParam(which = which, what = what, flag = flag147)
  problem_seq147 <- scanBam(file = bamPath, index = indexBam(bamPath), param = param)
  names(problem_seq147) <- NULL
  
  
  # FLAG 99
  
  param <- ScanBamParam(which = which, what = what, flag = flag99)
  problem_seq99 <- scanBam(file = bamPath, index = indexBam(bamPath), param = param)
  names(problem_seq99) <- NULL
  
  # FLAG 83 
  
  param <- ScanBamParam(which = which, what = what, flag = flag83)
  problem_seq83 <- scanBam(file = bamPath, index = indexBam(bamPath), param = param)
  names(problem_seq83) <- NULL
  #problem_seq83[[1]]$seq <- reverseComplement(problem_seq83[[1]]$seq)
  
  # FLAG 163
  
  param <- ScanBamParam(which = which, what = what, flag = flag163)
  problem_seq163 <- scanBam(file = bamPath, index = indexBam(bamPath), param = param)
  names(problem_seq163) <- NULL
  #problem_seq163[[1]]$seq <- reverseComplement(problem_seq163[[1]]$seq)
  
  # FLAG 16
  
  param <- ScanBamParam(which = which, what = what, flag = flag16)
  problem_seq16 <- scanBam(file = bamPath, index = indexBam(bamPath), param = param)
  names(problem_seq16) <- NULL
  #problem_seq163[[1]]$seq <- reverseComplement(problem_seq163[[1]]$seq)

  
  # Modify 83 AND 163
  
  
  problem_def <- list()
  problem_def$pos <- c(problem_seq147[[1]]$pos, problem_seq83[[1]]$pos, problem_seq99[[1]]$pos, problem_seq163[[1]]$pos, problem_seq16[[1]]$pos)
  problem_def$qwidth <- c(problem_seq147[[1]]$qwidth, problem_seq83[[1]]$qwidth, problem_seq99[[1]]$qwidth, problem_seq163[[1]]$qwidth, problem_seq16[[1]]$qwidth)
  problem_def$seq <- c(problem_seq147[[1]]$seq, problem_seq83[[1]]$seq, problem_seq99[[1]]$seq, problem_seq163[[1]]$seq, problem_seq16[[1]]$seq)
  
  return(problem_def)
  
  
}



