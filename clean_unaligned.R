
clean_unaligned <- function(align_data, binary_data) {

  col_names <- as.numeric(colnames(binary_data))
  binary_data2 <- binary_data
  for (i in 1:length(align_data)) {
  
    tmp_value <- align_data[[i]]@pattern@range
    binary_data2[i,] <- ifelse(start(tmp_value) <= col_names & col_names <= end(tmp_value), 'OK', 'NA' )
    binary_data[i,][which(binary_data2[i,] == 'NA')]  <- 2
    result <- binary_data
  }
  return(result)
}
