
together_gw <- function(df_end, df_exo, tss = NULL) {
  
  df_end$type <- rep('C', nrow(df_end))
  df_exo$type <- rep('G', nrow(df_exo))
  
  colnames(df_end) <- c('Position', 'Percentages', 'Legend')
  colnames(df_exo) <- c('Position', 'Percentages', 'Legend')
  
  
  df_merge <- rbind(df_end, df_exo)
  name_x <- "Reference sequence (pb)"
  if ( !is.null(tss)) {
    name_x <- "Distance to TSS (bp)"
    df_merge$Position <- df_merge$Position - as.integer(tss)
  }
  plot <- ggplot(data = df_merge, aes(x = Position, y = Percentages, group = Legend, colour = Legend)) + 
    geom_line(size=1.25) +
    geom_point(colour = "black", size = 1) +
    ylim(0, 100) +
    ylab("Average signal (%)") +
    xlab(name_x) +
    theme(plot.title = element_text(size=17)) +
    ggtitle("Overlap") +
    scale_colour_manual(values = c('black', 'red'),
                      labels = c('CpG methylation (%)', '1 - GpC methylation (%)'))
    
    # theme(legend.position="none")
  plot

}
