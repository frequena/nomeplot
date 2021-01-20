MethAlignNW_mod <- function (refSeq, QCdata, alignment) 
{
  seqName <- c()
  alignment <- c()
  files <- paste(QCdata$FILE)
  paths <- QCdata$PATH
  mat <- matrix(0, length(DNA_ALPHABET), length(DNA_ALPHABET))
  mat[1:4, 1:4] <- c(1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 
                     0, 0, 1)
  rownames(mat) <- colnames(mat) <- DNA_ALPHABET[1:length(DNA_ALPHABET)]
  lengthRef <- nchar(refSeq)
  cg_ref <- gregexpr("CG", refSeq)
  positionCGIRef <- sort(cg_ref[[1]][1:length(cg_ref[[1]])])
  methPos <- matrix(nrow = length(QCdata$PATH), ncol = length(cg_ref[[1]]))
  startEnd <- matrix(nrow = length(QCdata$PATH), ncol = 2)
  for (i in 1:length(paths)) {
    seqFileTemp <- paste(paths[i], "/", files[i], sep = "")
    seqTemp <- toupper(scan(seqFileTemp, what = "character", 
                            sep = "", quiet = TRUE))
    align <- pairwiseAlignment(refSeq, seqTemp, substitutionMatrix = mat, 
                               type = "local-global")
    startRef <- start(pattern(align))
    endRef <- end(pattern(align))
    startEnd_temp <- c(startRef, endRef)
    startEnd[i, ] <- startEnd_temp
    newRef <- paste(substr(refSeq, 1, startRef - 1), pattern(align), 
                    substr(refSeq, endRef + 1, lengthRef), sep = "")
    newSeq <- paste(paste(rep("-", startRef - 1), collapse = ""), 
                    subject(align), paste(rep("-", lengthRef - endRef)), 
                    sep = "")
    alignment[i] <- paste(subject(align))
    methPos_Temp <- cgMethFinder(newRef, newSeq)
    methPos[i, ] <- methPos_Temp
    seqName[i] <- files[i]
    cat("Alignment with ", seqName[i], " done", "\n")
  }
  if (missing(alignment) || alignment == FALSE) {
    mathylAlign <- list(seqName = seqName, methPos = methPos, 
                        positionCGIRef = positionCGIRef, startEnd = startEnd, 
                        lengthRef = nchar(refSeq))
  }
  else {
    mathylAlign <- list(seqName = seqName, alignment = alignment, 
                        methPos = methPos, positionCGIRef = positionCGIRef, 
                        startEnd = startEnd, lengthRef = nchar(refSeq))
  }
  return(mathylAlign)
}