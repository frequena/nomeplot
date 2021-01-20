window_df <- function(df, window_size = NULL) {
  if (is.null(window_size)) {
    size_window <- 120
  }
  else {
    size_window <- window_size
    
  }
  matrix_value <-  matrix(ncol = 3)
  matrix_value <- matrix_value[-1,]
  value_df <- as.integer(colnames(df))
  for (i in 2:(length(df)-1)) {
    diff_value_left <- value_df[i] - size_window/2
    diff_value_right <- value_df[i] + size_window/2
    number_left <- which.min(abs(value_df - diff_value_left))
    number_right <- which.min(abs(value_df - diff_value_right))
    if (diff_value_left > value_df[number_left] ) {
      number_left <- number_left + 1
    }
    if (diff_value_right < value_df[number_right]) {
      number_right <- number_right - 1
    }
    matrix_value <- rbind(matrix_value, c(value_df[number_left], value_df[number_right], 0))
  }
  number_right <- which.min(abs(value_df - (value_df[1] + (size_window/2))))
  number_left <- which.min(abs(value_df - (value_df[length(value_df)] -(size_window/2))))
  if ((value_df[length(value_df)] - size_window/2) > value_df[number_left]) {
    number_left <- number_left + 1
  }
  if ((value_df[1] + size_window/2) < value_df[number_right]) {
    number_right <- number_right - 1
  }

  matrix_value <- rbind(c(value_df[1], value_df[number_right], 0), matrix_value)
  matrix_value <- rbind(matrix_value, c(value_df[number_left], value_df[length(value_df)], 0))
  
  for (k in 1:(nrow(matrix_value))) {
    interval_value <- df[, which(colnames(df) == as.character(matrix_value[k,1])) : 
                           which(colnames(df) == as.character(matrix_value[k,2]))]
    mean_value <- mean(as.numeric(interval_value))
    matrix_value[k,3] <- mean_value
  }
  return(matrix_value)
}