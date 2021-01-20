
cleaning_df <- function(df) {
  
  df_loc <- df[c((nrow(df)-1), nrow(df)),]
  col_delete <- c()
  for (i in (1:(ncol(df_loc)-1))) {
    if (df_loc[1,][,i+1] == df_loc[1,][,i]) {
      col_delete <- append(col_delete, i)
    } 
  }
  df_loc <- subset(df_loc, select = -col_delete)
  col_delete <- c()
  for (i in 1:(ncol(df_loc)-1)) {
    if  (df_loc[2,][,i+1] == df_loc[2,][,i]) {
      col_delete <- append(col_delete, i+1)
    }
  }
  for (i in 1:(ncol(df_loc)-1)) {
      if  (df_loc[1,][,i] == df_loc[2,][,i]) {
          df_loc <- subset(df_loc, select = -i)
        }
     }
  df <- subset(df, select = as.numeric(names(df_loc)))
  return(df)
}