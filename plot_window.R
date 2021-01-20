plot_window <- function(df_exo, size_wd, name_seq, order_seq, loc_tss = NULL, type = NULL, names_column = NULL)
{
  n_row <- nrow(df_exo)
  df_exo <- 1 - df_exo
  # SCALE SIZE OF AXIS X LABELS
  size_legend <- 10
  if (n_row >= 20) {
    size_legend <- 8 
  }
  
  name_x <- "Reference sequence (pb)"
  
  if (!is.null(loc_tss)) {
    colnames(df_exo) <- (as.numeric(colnames(df_exo)) - as.numeric(loc_tss))
    name_x <- "Distance to TSS (bp)"
  }
  
  
  # Preparation data
  
  colnames(df_exo) <- colnames(names_column)
  df <- melt(df_exo)
  df$sequence <- name_seq[order_seq]
  df$sequence <- as.factor(df$sequence)
  
  df$sequence <- factor(df$sequence, levels= name_seq[order_seq] )
  # Preparation plot
  
  plot <- ggplot(df, aes(x = variable, y = value, color= sequence, group = sequence)) +
    geom_line() +
    xlab(name_x) +
    ylab(paste("Estimated nucleosome occupancy \n (1-average GpC methylation in", size_wd ,"bp regions")) +
    scale_fill_gradient(low="pink", high="blue") +
    theme(legend.title=element_blank(), 
          legend.text = element_text(size = size_legend),
          plot.title = element_text(size=16),
          axis.title=element_text(size=14,face="bold"),
          axis.text.x = element_text(angle = 90, hjust = 1))
  
  # Preparation of data
  
  df$variable <- as.factor(df$variable)
  df$sequence <- as.factor(df$sequence)
  # Generation of heatmap_plot
  
  heatmap_plot <- ggplot(df, aes(variable, sequence)) +
    geom_tile(aes(fill = value), color = "white") +
    scale_fill_gradient(name = paste("Estimated nucleosome occupancy \n (1-average GpC methylation \n in ", size_wd,  "bp regions"), 
                        low = "steelblue", high = "red") +
    ylab("NOMe-treated sequences") +
    xlab(name_x) +
    theme_bw() +
    theme( 
      legend.text = element_text(size = size_legend),
      plot.title = element_text(size=16),
      axis.title=element_text(size=14,face="bold"),
      axis.text.x = element_text(angle = 90, hjust = 1))
  if (type == "plot") {
    print(plot)
  } else {
    print(heatmap_plot)
  }
}
