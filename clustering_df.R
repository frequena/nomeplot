
clustering_df <- function(df_translate, r_problem, names_cond = NULL) {
  if (!is.null(names_cond)) {
  df_cluster_order <- hclust(dist(df_translate))
  rownames(df_translate) <- names(r_problem)
  df_cluster_dend <- hclust(dist(df_translate))
  result <- list("dend" = df_cluster_dend, "order" = df_cluster_order)
  } else {
    result <-df_cluster_dend <- hclust(dist(df_translate))
  }
  return(result)
}