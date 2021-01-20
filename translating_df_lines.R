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

