read_seq_pro <- function(m_sequence) {
  data_read <- toupper(readLines(m_sequence))
  name_seq <- data_read[1]
  data_read <- data_read[-1]
  result <- list("SEQ" = DNAString(paste(data_read, sep = "\n", collapse = "")), "NAME" = name_seq)
return(result)
}