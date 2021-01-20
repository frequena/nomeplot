read_seq_ref <- function(m_sequence) {
  data_read <- toupper(readLines(m_sequence))
  name_seq <- data_read[1]
  data_read <- data_read[-1]
  result <- DNAString(paste(data_read, sep = "\n", collapse = ""))
return(result)
}