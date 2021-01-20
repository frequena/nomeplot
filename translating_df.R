
translating_df <- function(df_exo, size_wd = NULL, pos_nucleosome = TRUE) {
  number_row <- nrow(window_df(df_exo[1,], window_size = size_wd))
  df_translate <- data.frame(matrix(nrow = nrow(df_exo), ncol = number_row))
  for (i in 1:nrow(df_exo)) {
    Temp_value <- window_df(df_exo[i,], size_wd)
    df_translate[i,] <- Temp_value[,3] 
  }
  colnames(df_translate) <- c(1:length(df_translate))
  if (pos_nucleosome == TRUE) {
    df_translate <- 1 - df_translate
  }
  
  return(df_translate)
}


translating_df_lines <- function(df_exo, df_exo2, size_wd = NULL) {
  colnames(df_exo) <- colnames(df_exo2)
  df_translate_lines <- data.frame(matrix(nrow = (nrow(df_exo)+2), ncol = ncol(df_exo)),
                                   check.names = FALSE)
  df_translate_lines[nrow(df_exo)+1,] <- window_df(df_exo[1,], size_wd)[,1]
  df_translate_lines[nrow(df_exo)+2,] <- window_df(df_exo[1,], size_wd)[,2]
  for (i in 1:nrow(df_exo)) {
    df_translate_lines[i,] <- window_df(df_exo[i,], size_wd)[,3]
  }
  df_translate_lines <- cleaning_df(df_translate_lines)
  return(df_translate_lines)
}


translating_df_plot <- function(df_exo, size_wd = NULL) {
  df_translate_plot <- data.frame(matrix(nrow = nrow(df_exo), ncol = ncol(df_exo)),
                                  check.names = FALSE)
  df_exo <- 1 - df_exo
  for (i in 1:nrow(df_exo)) {
    df_translate_plot[i,] <- window_df(df_exo[i,], size_wd)[,3]
  }
  return(df_translate_plot)
}

