nucleoplot_plot <- function(x, tss = NULL, unmeth = NULL, n_occupancy = NULL, box = NULL ) {
  
  require(ggplot2)
  n_result <- c()
  n_row <- nrow(x)
  n_col <- ncol(x)
  col_names <- as.numeric(colnames(x))
  value_meth <- 1
  if (unmeth == TRUE) {
    value_meth <- 0
  }
  for (i in 1:n_col) {
    n_meth <- length(as.numeric(x[,i][x[,i] == value_meth]))
    n_result <- append(n_result,n_meth)
  }
  n_result <- (n_result / n_row) * 100
  n_result <- format(round(n_result, 2), nsmall = 2)
  name_x <- as.numeric(colnames(x))
  if (!is.null(tss)) {
    name_x <- as.numeric(colnames(x)) - tss
  }
  
  df_analyses_exo <- data.frame("Names" =name_x,
                                "Values" = as.numeric(n_result))
  if ( unmeth == TRUE && n_occupancy == TRUE ) {
    name_x <- "Nucleosome occupancy \n (% of unmethylated GpC sites)"
    title_mod <- "Nucleosome occupancy \n (% of unmethylated GpC)"
    color_mod <- "red"
  }
  # REVISAR!!
  if ( unmeth == TRUE && n_occupancy == FALSE) {
    name_x <- "% of methylated GpC sites"
  }
  if ( unmeth == FALSE) {
    title_mod <- "Endogenous DNA methylation \n (% of methylated CpG)"
    name_x <- "% of methylated CpG sites"
    color_mod <- "black"
  }
  name_y <- "Reference sequence (pb) "
  
  if ( !is.null(tss)) {
    name_y <- "Distance to TSS (bp)"
  } 
  
  if (is.null(box) == TRUE ) {
    pd <- position_dodge(.1)
    plot <- ggplot(data = df_analyses_exo, aes(x = Names, y = Values, group = 1)) + 
      geom_line(position=pd, colour = color_mod, size=1.25) +
      geom_point(size = 1) +
      ylim(0, 100) +
      ylab(name_x) +
      xlab(name_y) +
      ggtitle(title_mod)+
      theme(plot.title = element_text(size=17))
  } else {
    
    pd <- position_dodge(.1)
    plot <- ggplot(data = df_analyses_exo, aes(x = Names, y = Values, group = 1)) + 
      geom_line(position=pd, colour = color_mod, size=1.25) +
      geom_point(size = 1) +
      ylim(0, 100) +
      ylab(name_x) +
      xlab(name_y) +
      ggtitle(title_mod)+
      theme(plot.title = element_text(size=17))
    
    
    
  }
  return(plot)
}
