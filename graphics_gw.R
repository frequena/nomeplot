

graphics_gw <- function(data_info, cpg = NULL,  table = NULL) {
  
  
    colnames(data_info) <- as.character(as.numeric(colnames(data_info)) - 2)
    
    
  
    discard_dupl <- c()
    
  for (i in 1:ncol(data_info)) {
    
     if (sum(!duplicated(data_info[,i]))==1 & data_info[,i][1] == '2') {
       
       discard_dupl <- append(discard_dupl, as.numeric(i))
       
    
  }
  }
    if (length(discard_dupl > 0)) {
      
      data_info <- data_info[,-discard_dupl]
    }
    
   
  lv_perc <- c()
 
  for (i in 1:ncol(data_info)) {
    
    temp_data <- as.data.frame(table(data_info[,i]))
    
  
    temp_zero <- as.numeric(filter(temp_data, Var1 == 0)[,2])
    if (length(temp_zero) == 0) {
      temp_zero <- length(temp_zero)
    }
   

    temp_one <- as.numeric(filter(temp_data, Var1 == 1)[,2])
    if (length(temp_one) == 0) {
      temp_one <- length(temp_one)
    }
     
  
    temp_value <- temp_one /(temp_zero + temp_one)
    lv_perc <- round(append(lv_perc, temp_value), 2)
  }
  data_perc <- data.frame('names' = as.numeric(colnames(data_info)), 'values' = lv_perc * 100 )
  
  
  # Plot title
  
  name_x <- "Nucleosome occupancy \n (% of unmethylated GpC sites)"
  name_y <- "Reference sequence (pb) "
  title_mod <- "Nucleosome occupancy \n (% of unmethylated GpC)"
  color_mod <- "red"
  
  if (cpg == TRUE) {
    title_mod <- "Endogenous DNA methylation \n (% of methylated CpG)"
    name_y <- "Reference sequence (pb) "
    name_x <- "% of methylated CpG sites"
    # color_mod <- "brown2"
    color_mod <- 'black'
    
    
  }
  if (cpg == FALSE) {
    data_perc$values <- 100 - data_perc$values
  }
  pd <- position_dodge(.1)
  plot <- ggplot(data = data_perc, aes(x = names, y = values, group = 1)) + 
    geom_line(position=pd, colour = color_mod, size=1.25) +
    geom_point(size = 1) +
    ylim(0, 100) +
    ylab(name_x) +
    xlab(name_y) +
    ggtitle(title_mod)+
    theme(plot.title = element_text(size=17))
  
  plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  if (!is.null(table)) {
    
    colnames(data_perc) <- c("GpC.sites", '% of unmethylated GpC sites')
    if (cpg == TRUE) {
      
      colnames(data_perc) <- c("CpG.sites", "% of methylated CpG")
    }
    
    plot <- data_perc
    
  }
  
  return(plot)
  
}
