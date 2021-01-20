
nomeplot_merge <- function(resume_table, tss = NULL) {
  resume_table_end <- resume_table$END
  resume_table_end$Legend <- "CpG methylation"
  resume_table_exo <- resume_table$EXO[,c(1,3)]
  resume_table_exo$Legend <- "Nucleosome occupancy"
  names(resume_table_exo)[1:2] <- c("Position", "Percentages")
  names(resume_table_end)[1:2] <- c("Position", "Percentages")
  
  df_merge <- rbind(resume_table_end, resume_table_exo)
  df_merge[,2] <- df_merge[,2] * 100
  name_x <- "Reference sequence (pb)"
  if ( !is.null(tss)) {
    name_x <- "Distance to TSS (bp)"
    df_merge$Position <- df_merge$Position - as.integer(tss)
  }
  plot <- ggplot(data = df_merge, aes(x = Position, y = Percentages, group = Legend, colour = Legend)) + 
  geom_line(size=1.25) +
  geom_point(colour = "black", size = 1) +
  ylim(0, 100) +
  ylab("1-GpC methylation / CpG methylation") +
  xlab(name_x) +
  scale_color_manual(values=c("black", "red")) +
  theme(legend.title=element_blank() , legend.text = element_text(size = 14), 
        plot.title = element_text(size=17)) +
  ggtitle("Overlap")
return(plot)
}
