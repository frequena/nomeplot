
read_master <- function(m_sequence) {
  require(Biostrings)
  if ( m_sequence != "") {
    data_read <- toupper(readLines(m_sequence))
    data_name <- data_read[[1]]
    m_seq <- DNAString(data_read[[2]])
  }
  else {
    cat("Not master sequence detected. Program stopped.")
  }
  return(m_seq)
}