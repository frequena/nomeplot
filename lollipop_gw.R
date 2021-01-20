lollipop_gw <- function(df_exo, tss = NULL, col_circles = NULL, original_label = NULL, gpc = NULL) {
  # set.seed(123)
  # df_exo <- gpc_gw_def
  # gpc <- F
  # tss <- NULL
  # require(grid)
  
  # if (gpc == T) {
  #   
  #   colnames(df_exo) <- as.character(as.numeric(colnames(df_exo))
  # }
  
  #### DISCARD ROWS WITH NA (NO FUNCIONA CUANDO APLICO EL QUITAR LOS 2 FULL ROWS)
  discard_rows_na <- c()
  
  for (i in 1:nrow(df_exo)) {
    
   temp_na <- grep('NA', df_exo[i,])
    if (length(temp_na) > 0) {
      discard_rows_na <- append(discard_rows_na, i)
    }
  }
  if (length(discard_rows_na) > 0) {
  df_exo <- df_exo[-discard_rows_na,]
  }
  
  # ADJUSTMENT OF DIMENSION
  
  # No sé que es esto la verdad...
  if (!is.null(original_label)) {
    
    colnames(df_exo) <- as.numeric(colnames(df_exo)) + original_label
  }
  
  # if (nrow(df_exo) > 46) {
  # 
  #   #df_exo <- df_exo[1:46,]
  #   df_exo <- df_exo[sample(1:nrow(df_exo), 46),]
  #   df_exo <- df_exo[ order(as.numeric(row.names(df_exo))), ]
  # }
  # 
  # ncol_problem <- nrow(df_end)
  
  rownames(df_exo) <- as.character(1:nrow(df_exo))
  # Order rows by coverage
  if (nrow(df_exo) > 46) {
  number_two <- c()
  for (i in 1:nrow(df_exo)) {
    
    
    number_two <- append(number_two, rowSums(df_exo[i,] == 2))
    
  }

  number_two <- sort(number_two, decreasing = F)
  df_exo <- df_exo[as.numeric(names(number_two)),]
  if (nrow(df_exo) < 45) {
    number_choosen <- nrow(df_exo)
  } else {
    number_choosen <- 45
  }
  df_exo <- df_exo[sample(1:nrow(df_exo),number_choosen),] # Random selection
  # df_exo <- df_exo[1:45,] # Selection
  }
  
  ##
  # size_circles <- 0.008
  adjust_y_end <- 1
  adjust_y_exo <- 1
  
  # if (ncol_problem <= 20) {
     size_circles <- 0.008
     adjust_y_exo <- 0.019 
  #   adjust_y_end <- 0.019
  #   height_line <- 0.015
  # }   else if (ncol_problem >= 21 & ncol_problem <= 28) {
  #   size_circles <- 0.007
  #   adjust_y_exo <- 60 
  #   adjust_y_end <- 0.015
  #   height_line <- 0.013
  #   
  # }   else if (ncol_problem >= 29 & ncol_problem <= 36) {
  #   size_circles <- 0.004
  #   adjust_y_exo <- 77 
  #   adjust_y_end <- 0.012
  #   height_line <- 0.011
  #   
  # }   else if (ncol_problem >= 37 & ncol_problem <= 44) {
  #   size_circles <- 0.003
  #   adjust_y_exo <- 90 
  #   adjust_y_end <- 0.009
  #   height_line <- 0.009
  #   
  # }   else {
  #   size_circles <- 0.002
  #   adjust_y_exo <- 85  
  #   adjust_y_end <- 0.010
  #   height_line <- 0.007
  #   
  # }
  
  # COLOR_CIRCLES
  
  if (is.null(col_circles)) {
    color_circles <- "blue"
  }
  else {
    color_circles <- col_circles
  }  
  
  # COLOR_LINES
  dict_color <- list("red" = rgb(1,0,0,0.8) , 
                     "blue" =rgb(0,0,1,0.8) ,
                     "black" =  rgb(0,0,0,0.8), 
                     "cyan" =  rgb(0,1,1,0.8),
                     "magenta" =  rgb(1,0,1,0.8),
                     "green" =  rgb(0,1,0,0.8))
  
 # General values
  
  n_row <- nrow(df_exo)
  n_col <- ncol(df_exo)
  
  # Change of order
  
  # df_exo <- df_exo[order(-1:-n_row),]
  # rownames(df_exo) <- c(1:n_row)
  # df_translate <- df_translate[order(-1:-n_row),]
  # rownames(df_translate) <- c(1:nrow(df_translate))
  
  # Generation of vector with positions of every GpC site
  
  location_end <- as.numeric(colnames(df_exo))
  #location_end <- location_end - interval_min
  
  # Generation of every position in relation with TSS
  
  location_mod <- ((location_end * 0.95)/max(location_end))
  
  # CIRCLES #
  
  #  ((j*0.95)/adjust_y_exo)
  # (0.885 - (j-1)*adjust_y_end)
  
    for (j in n_row:1){
      for (i in 1:n_col){
        k = location_mod[i]
        value_translate = as.numeric(df_exo[j,][,i])
        if (value_translate == 1){
          result  <- grid.circle(x=k, y=(0.885 - (j-1)*adjust_y_exo), r=size_circles, name="circles",  
                                 gp=gpar(fill=color_circles))
        }
        if (value_translate == 0){
          result  <- grid.circle(x=k, y=(0.885 - (j-1)*adjust_y_exo), r=size_circles,                          
                                 name="circles")
        }
      }
      
    }
 
  for (l in 1:n_col) {
   # if (is.null(tss)) {
      # result <- grid.text((location_end[l]), x = location_mod[l], y = 0.925, 
      #                     gp=gpar(fontfamily="Helvetica", fontsize = 9), name="font",
      #                     rot = 90,
      #                     just = c("centre","centre"))
      
      # result <- grid.text((colnames(df_exo)[l]), x = location_mod[l], y = 0.925, 
      #                     gp=gpar(fontfamily="Helvetica", fontsize = 7), name="font",
      #                     rot = 90,
      #                     just = c("centre","centre"))
    
      location_end_tss <- location_end - tss
      result <- grid.text(location_end_tss[l], x = location_mod[l], y = 0.925, 
                          gp=gpar(fontfamily="Helvetica", fontsize = 9), name="font", rot = 90,
                          just = c("centre","centre"))
    
    }

  # # DRAWN OF ARROW TSS #
  #     tss_x <- ((tss)*0.95)/result_max
  #     result <- grid.lines(x = unit(c(tss_x, tss_x), "npc"),
  #                          y = unit(c(0.94, 0.98), "npc"),
  #                          gp=gpar(col= 1, lwd = 4))
  #     result <- grid.lines(x = unit(c(tss_x, (tss_x + 0.065)), "npc"),
  #                          y = unit(c(0.98, 0.98), "npc"),
  #                          arrow = arrow(angle = 30, length = unit(0.15, "inches"),
  #                                        ends = "last", type = "open"),
  #                          gp=gpar(col= 1, lwd = 4))
  #  
 
}


