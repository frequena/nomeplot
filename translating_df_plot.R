translating_df_plot <- function(df_exo, size_wd = NULL) {
  df_translate_plot <- data.frame(matrix(nrow = nrow(df_exo), ncol = ncol(df_exo)),
                                   check.names = FALSE)
  df_exo <- 1 - df_exo
  for (i in 1:nrow(df_exo)) {
    df_translate_plot[i,] <- window_df(df_exo[i,], size_wd)[,3]
  }
  return(df_translate_plot)
}

