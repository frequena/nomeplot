lollipop_svg_doble_lines4 <- function(df_end, df_exo, df_translate, nucleosome_mean, tss = NULL
                                      , lines = NULL, window_size = NULL,
                                      col_lines = NULL, col_circles = NULL) {
  require(grid)
  
  df_exo <- df_exo %>% na.omit()
  df_end <- df_end %>% na.omit()
  
  
  # NUMBER OF DRAWN LINES
  
  number_lines <- round(as.numeric(colnames(df_exo)[length(df_exo)]) / 200)
  
  # ADJUSTMENT OF DIMENSION
  
  ncol_problem <- nrow(df_end)
  size_circles <- 1
  adjust_y_end <- 1
  adjust_y_exo <- 1

  if (ncol_problem <= 20) {
    size_circles <- 0.009
    adjust_y_exo <- 50 
    adjust_y_end <- 0.019
    height_line <- 0.015
  }   else if (ncol_problem >= 21 & ncol_problem <= 28) {
    size_circles <- 0.007
    adjust_y_exo <- 60 
    adjust_y_end <- 0.015
    height_line <- 0.013
    
  }   else if (ncol_problem >= 29 & ncol_problem <= 36) {
    size_circles <- 0.004
    adjust_y_exo <- 77 
    adjust_y_end <- 0.012
    height_line <- 0.011
    
  }   else if (ncol_problem >= 37 & ncol_problem <= 44) {
    size_circles <- 0.003
    adjust_y_exo <- 90 
    adjust_y_end <- 0.009
    height_line <- 0.009
    
  }   else {
    size_circles <- 0.002
    adjust_y_exo <- 85  
    adjust_y_end <- 0.010
    height_line <- 0.007
    
  }
  
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
  
  if (is.null(col_lines)) {
    color_lines <- rgb(1,0,0,0.8)
  } else {
    color_lines <- dict_color[col_lines][[1]]
  }
  
  # MAXIMUM VALUE OF DATA.FRAME
  
  max_end <- max(as.numeric(colnames(df_end)))
  max_exo <- max(as.numeric(colnames(df_exo)))
  
  result_max <- max_exo
  if (max_end >= max_exo) {
    result_max <- max_end
  }
  
  
  # DF_EXO
  
  # General values
  
  n_row <- nrow(df_exo)
  n_col <- ncol(df_exo)
  
  # Change of order
  
  df_exo <- df_exo[order(-1:-n_row),]
  rownames(df_exo) <- c(1:n_row)
  df_translate <- df_translate[order(-1:-n_row),]
  rownames(df_translate) <- c(1:nrow(df_translate))
  
  # Generation of vector with positions of every GpC site
  
  location_end <- as.numeric(colnames(df_exo))
  
  # Generation of every position in relation with TSS
  
  location_mod <- ((location_end * 0.95)/result_max)
  
  # LINEAS NULL (SOLO CIRCLES) #
  
  if (is.null(lines)) {
    for (j in n_row:1){
      for (i in 1:n_col){
        k = location_mod[i]
        value_translate = as.numeric(df_exo[j,][,i])
        value_translate <<- value_translate
        if (value_translate == 1){
          result  <- grid.circle(x=k, y=((j*0.95)/adjust_y_exo), r=size_circles, name="circles",  
                                 gp=gpar(fill=color_circles))
        }
        if (value_translate == 0){
          result  <- grid.circle(x=k, y=((j*0.95)/adjust_y_exo), r=size_circles,                          
                                 name="circles")
        }
      } 
    }
  }
  else {
    
    # DRAWN OF LINES
    
    # 1. Generation of df_translate of drawn lines
    
    df_translate_lines <<- translating_df_lines(df_translate, df_exo , window_size)
    
    # 2. Configuration of total distance of work
    
    total_distance <- max(location_mod) - min(location_mod)
    
    # 3. Configuration of relation npc - bp
    
    relation_npc_nt <- total_distance/max(location_end)
    
    # 4. Start of every line an drawn
    
    count_lines <- 0
    
    for (j in n_row:1){
      for (i in 1:n_col){
        k = location_mod[i]
        value_translate = as.numeric(df_exo[j,][,i])
        if (value_translate == 1){
          result  <- grid.circle(x=k, y=((j*0.95)/adjust_y_exo), r=size_circles, name="circles",  
                                 gp=gpar(fill=color_circles))
        }
        if (value_translate == 0){
          result  <- grid.circle(x=k, y=((j*0.95)/adjust_y_exo), r=size_circles,                          
                                 name="circles")
        }
      }
    }
    for (j in n_row:1){
      count_lines <- 0
      last_position <- 0
      for (i in 1:(ncol(df_translate_lines)-1)) {
        if (count_lines >= number_lines) break
        if (df_translate_lines[j,][,i] >= nucleosome_mean & 
               last_position <= df_translate_lines[(nrow(df_translate_lines)-1),][,i]) {
          
          last_position <- df_translate_lines[(nrow(df_translate_lines)),][,i]
          
          # Length of line
          length_line <- ((df_translate_lines[(nrow(df_translate_lines)-1),][,i] * 0.95)/result_max) -
            ((df_translate_lines[(nrow(df_translate_lines)),][,i] * 0.95)/result_max)
          # Start of line
          
          start_line <- ((df_translate_lines[(nrow(df_translate_lines)-1),][,i] * 0.95)/result_max) + 0.02
          if ((length_line + start_line) >= location_mod[length(location_mod)]) {
            length_line <- length_line - ((start_line+length_line)- location_mod[length(location_mod)])
          }

          grid.rect(x=(start_line), y = (((j*0.95)/adjust_y_exo)),
                    width=unit(length_line, "npc"),
                    height=unit(height_line, "npc"),
                    gp = gpar(fill = color_lines , col = color_lines),
                    just = "right")
          count_lines <- count_lines + 1
        }
      }
    } 
  }
  # ELABORATION GpC LABEL #
  result <- grid.text("GpC", x = 0.01, y = (((n_row/2)*0.95)/adjust_y_exo) , 
                      gp=gpar(fontfamily="Helvetica", fontsize = 15), name="font", rot = 90,
                      just = c("centre","centre"))
  
  for (l in 1:n_col) {
    if (is.null(tss)) {
      result <- grid.text((location_end[l]), x = location_mod[l], y = (((n_row+2)*0.95)/adjust_y_exo), 
                          gp=gpar(fontfamily="Helvetica", fontsize = 9), name="font",
                          rot = 90, just = c("centre", "centre"))
    }
    else {
      location_end_tss <- location_end - tss
      result <- grid.text(location_end_tss[l], x = location_mod[l], y = (((n_row+2)*0.95)/adjust_y_exo), 
                          gp=gpar(fontfamily="Helvetica", fontsize = 9), name="font", rot = 90,
                          just = c("centre", "centre"))
    }
  }
  
  # DF_END
  n_row <- nrow(df_end)
  n_col <- ncol(df_end)
  location_end <- as.numeric(colnames(df_end))
  location_end_tss <- location_end - tss
  location_mod <- (location_end*0.95)/result_max
  for (j in 1:n_row){
    for (i in 1:n_col){
      k = location_mod[i]
      value_translate = as.numeric(df_end[j,][,i])
      if (value_translate == 1){
        result  <- grid.circle(x=k, y=(0.885 - (j-1)*adjust_y_end), r=size_circles, name="circles",  
                               gp=gpar(fill="black"))
      }
      if (value_translate == 0){
        result  <- grid.circle(x=k, y=(0.885 - (j-1)*adjust_y_end), r=size_circles,                          
                               name="circles")
      }
    }
  }
  for (l in 1:n_col) {
    if (is.null(tss)) {
      result <- grid.text((location_end[l]), x = location_mod[l], y = 0.925, 
                          gp=gpar(fontfamily="Helvetica", fontsize = 9), name="font",
                          rot = 90,
                          just = c("centre","centre"))
    }
    else {
      location_end_tss <- location_end - tss
      result <- grid.text(location_end_tss[l], x = location_mod[l], y = 0.925, 
                          gp=gpar(fontfamily="Helvetica", fontsize = 9), name="font", rot = 90,
                          just = c("centre","centre"))
      # DRAWN OF ARROW TSS #
      tss_x <- ((tss)*0.95)/result_max
      result <- grid.lines(x = unit(c(tss_x, tss_x), "npc"),
                           y = unit(c(0.94, 0.98), "npc"),
                           gp=gpar(col= 1, lwd = 4))
      result <- grid.lines(x = unit(c(tss_x, (tss_x + 0.065)), "npc"),
                           y = unit(c(0.98, 0.98), "npc"),
                           arrow = arrow(angle = 30, length = unit(0.15, "inches"),
                                         ends = "last", type = "open"),
                           gp=gpar(col= 1, lwd = 4))
    }
  }
  # ELABORATION OF CpG LABEL #
  result <- grid.text("CpG", x = 0.01, y = (0.925 - ((n_row/2)-1)*adjust_y_end) , 
                      gp=gpar(fontfamily="Helvetica", fontsize = 15), name="font", rot = 90,
                      just = c("centre","centre"))
}
