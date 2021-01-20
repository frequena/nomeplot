lollipop_svg_doble_lines2 <- function(df_end, df_exo, df_translate, nucleosome_mean, tss = NULL
                                     , lines = NULL, col_lines = NULL, col_circles = NULL) {
  require(grid)
  
  # COLOR_CIRCLES
  
  if (is.null(col_circles)) {
    color_circles <- "blue"
  }
  else {
    color_circles <- col_circles
  }  
  
  # CREATION LIST (COLOR - RGB)
  dict_color <- list("red" = rgb(1,0,0,0.8) , 
                     "blue" =rgb(0,0,1,0.8) ,
                     "black" =  rgb(0,0,0,0.8), 
                     "cyan" =  rgb(0,1,1,0.8),
                     "magenta" =  rgb(1,0,1,0.8),
                     "green" =  rgb(0,1,0,0.8))
  
  # COLOR_LINES
  
  if (is.null(col_lines)) {
    color_lines <- rgb(1,0,0,0.8)
  } else {
    color_lines <- dict_color[col_lines][[1]]
  }
  
  # DF_EXO
  n_row <- nrow(df_exo)
  n_col <- ncol(df_exo)
  df_exo <- df_exo[order(-1:-n_row),]
  rownames(df_exo) <- c(1:n_row)
  df_translate <- df_translate[order(-1:-n_row),]
  rownames(df_translate) <- c(1:nrow(df_translate))
  location_end <- as.numeric(colnames(df_exo))
  location_end_tss <- location_end - tss
  max_position <- max(location_end)
  location_mod <- (location_end*0.95)/max_position
  
  # HIPOTESIS #
  total_distance <- max(location_mod) - min(location_mod)
  vector_column <- as.numeric(colnames(df_exo))
  random_value <- sort(c(vector_column[1], sample(vector_column[-c(1, length(vector_column))],
                                          ncol(df_translate)),vector_column[length(vector_column)]))
  dist_values <- diff(random_value)
  dist_values <- dist_values[-length(dist_values)]
  correlation_value <- total_distance / sum(dist_values)
    
  # LINEAS #
  length_line <- dist_values * correlation_value
  start_line <- c(min(location_mod)) - 0.005
  for (i in 1:(ncol(df_translate)-1)) {
    sum_distance <- c(min(location_mod) - 0.005) + sum(length_line[1:i])
    start_line <- as.vector(append(start_line, sum_distance))
  }
  
  # ADJUSTED LIMITS #
  length_line[length(length_line)-1] <- length_line[length(length_line)-1] +0.003
  start_line[length(start_line)] <- start_line[length(start_line)] + 0.003
  
  # LINEAS #
  
  if (is.null(lines)) {
    for (j in n_row:1){
      for (i in 1:n_col){
        k = location_mod[i]
        value_translate = as.numeric(df_exo[j,][,i])
        if (value_translate == 1){
          result  <- grid.circle(x=k, y=((j*0.95)/50), r=0.01, name="circles",  
                                 gp=gpar(fill=color_circles))
        }
        if (value_translate == 0){
          result  <- grid.circle(x=k, y=((j*0.95)/50), r=0.01,                          
                                 name="circles")
        }
      } 
    }
  }
  else {
    for (j in n_row:1){
      for (i in 1:n_col){
        k = location_mod[i]
        value_translate = as.numeric(df_exo[j,][,i])
        if (value_translate == 1){
          result  <- grid.circle(x=k, y=((j*0.95)/50), r=0.01, name="circles",  
                                 gp=gpar(fill=color_circles))
        }
        if (value_translate == 0){
          result  <- grid.circle(x=k, y=((j*0.95)/50), r=0.01,                          
                                 name="circles")
        }
      }
      for (i in 1:ncol(df_translate)) {
        if (df_translate[j,][,i] >= nucleosome_mean) {
          grid.rect(x=(start_line[i]), y = (((j*0.95)/50)),
                    width=unit(length_line[i], "npc"),
                    height=unit(0.015, "npc"),
                    gp = gpar(fill = color_lines , col = color_lines),
                    just = "left") 
        }
      }
    } 
  }
  
  for (l in 1:n_col) {
    if (is.null(tss)) {
      result <- grid.text(location_end[l], x = location_mod[l], y = (((n_row+2)*0.95)/50), 
                          gp=gpar(fontfamily="Helvetica", fontsize = 9), name="font")
    }
    else {
      location_end_tss <- location_end - tss
      result <- grid.text(location_end_tss[l], x = location_mod[l], y = (((n_row+2)*0.95)/50), 
                          gp=gpar(fontfamily="Helvetica", fontsize = 9), name="font", rot = 90,
                          just = c("centre", "centre"))
    }
  }
  # DF_END
  n_row <- nrow(df_end)
  n_col <- ncol(df_end)
  location_end <- as.numeric(colnames(df_end))
  location_end_tss <- location_end - tss
  max_position <- max(as.numeric(colnames(df_exo)))
  location_mod <- (location_end*0.95)/max_position
  for (j in 1:n_row){
    for (i in 1:n_col){
      k = location_mod[i]
      value_translate = as.numeric(df_end[j,][,i])
      if (value_translate == 1){
        result  <- grid.circle(x=k, y=(0.925- (j-1)*0.019), r=0.01, name="circles",  
                               gp=gpar(fill="black"))
      }
      if (value_translate == 0){
        result  <- grid.circle(x=k, y=(0.925- (j-1)*0.019), r=0.01,                          
                               name="circles")
      }
    }
  }
  for (l in 1:n_col) {
    if (is.null(tss)) {
      result <- grid.text(location_end[l], x = location_mod[l], y = 0.965, 
                          gp=gpar(fontfamily="Helvetica", fontsize = 9), name="font")
    }
    else {
      location_end_tss <- location_end - tss
      result <- grid.text(location_end_tss[l], x = location_mod[l], y = 0.965, 
                          gp=gpar(fontfamily="Helvetica", fontsize = 9), name="font", rot = 90,
                          just = c("centre","centre"))
    }
  }
}