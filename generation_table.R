generation_table <- function(dataframe_end, dataframe_exo, tss = NULL) {
  
  cpg_position <- as.integer(colnames(dataframe_end))
  gpc_position <-  as.integer(colnames(dataframe_exo))
  methy_cpg <-  round(as.numeric(unname(colSums(dataframe_end))/nrow(dataframe_end)), 2)
  methy_gpc <-  round(as.numeric(unname(colSums(dataframe_exo))/nrow(dataframe_exo)), 2)
  unmethy_gpc <- 1 - methy_gpc
  if ( !is.null(tss)) {
    cpg_position <- cpg_position - as.numeric(tss)
    gpc_position <- gpc_position - as.numeric(tss)
  }
  dataframe_cpg <- data.frame("CpG position" = cpg_position, "% of methylated CpG" = methy_cpg,
                              check.names = FALSE, row.names = NULL)
  dataframe_gpc <- data.frame("GpC position" = gpc_position, "% of methylated GpC" = methy_gpc,
                              "Nucleosome occupancy (% of unmethylated GpC)" = unmethy_gpc,
                              check.names = FALSE, row.names = NULL)
  rownames(dataframe_cpg) <- c()
  row.names(dataframe_gpc) <- c()

  result <- list("END" = dataframe_cpg, "EXO" = dataframe_gpc )
  return(result)
}