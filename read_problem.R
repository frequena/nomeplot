

read_problem <- function(p_sequences) {
  require(Biostrings)
  require(seqinr)
  data_read <- read.fasta(p_sequences, as.string = TRUE)
  info <- data.frame(lapply(data_read, nchar))
  info_data <- data.frame("Names problem sequences" = colnames(info), "Length(pb)" = as.numeric(info[1,]))
  seq_set <- c()
  seq <- ""
  for (i in 1:length(data_read)) {
    seq <- ""
    seq_set <- append(data_read[[i]][1], seq_set)
  }
  dna_set <- rev(DNAStringSet(seq_set))
  names_sequences <<- colnames(info)
  dna_set@ranges@NAMES <- names_sequences
  return(dna_set)
  print(names_sequences)
}
