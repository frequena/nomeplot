shinyServer(function(input, output, session) {
  
  source("read_problem.R")
  source("read_master.R")
  source("qc_seq.R")
  source("del_gcg_end.R")
  source("del_gcg_exo.R")
  source("align_seq.R")
  source("MethAlignNW_exo.R")
  source("MethAlignNW_mod.R")
  source("cgMethFinder_exo.R")
  source("cgMethFinder_mod.R")
  source("translating_df.R")
  source("window_df.R")
  source("clustering_df.R")
  source("nucleoplot_plot.R")
  source("primer_seq.R")
  source("read_seq_ref.R")
  source("read_seq_pro.R")
  source("lollipop_svg_lines.R")
  source("lollipop_svg_doble_lines4.R")
  source("generation_table.R")
  source("nomeplot_merge.R")
  source("multiplot.R")
  source("plot_window.R")
  source("plot_window_bi.R")
  source("cleaning_df.R")
  source("lollipop_bi.R")
  source('lollipop_gw.R')
  source('graphics_gw.R')
  source('collect_seq.R')
  source('vcf_discard.R')
  source('plotviz.R')
  source('getting_ref.R')
  source('together_gw.R')

  library(Biostrings)
  library(rmarkdown)
  library(reshape2)
  library(seqinr)
  library(BiocGenerics)
  library(Gviz)
  library(GenomicRanges)
  library(Rsamtools)
  library(tidyverse)
  library(shinydashboard)
  library(shinydashboardPlus)
  library(shinyWidgets)
  library(prettydoc)
  library(cowplot)
  library(rsconnect)
  library(PKI)
  library(RCurl)
  library(gridBase)
  library(gridExtra)
  library(gridSVG)
  library(grid)
  library(import)
  library(dplyr)
  library(shinycssloaders)
  library(httr)
  library(stringr)
  library(shinyjs)
  library(magrittr)

  options(showHeadLines = 100)
  
  
  
  observeEvent(input$tab == 'nome', {
    
    #
    shinyjs::reset("file_reference_bi")
    shinyjs::reset("files_problem_bi")
    shinyjs::reset("fw_primer_bi")
    shinyjs::reset("rv_primer_bi")
    shinyjs::reset("primers_bi")
    shinyjs::reset("num_tss_bi")
    shinyjs::reset("tss_bi")
    
    
    #
    shinyjs::reset('file_reference_gw')
    shinyjs::reset('chr_select')
    shinyjs::reset('start_interval')
    shinyjs::reset('end_interval')
    shinyjs::reset('run')
    shinyjs::reset('bam_gw')
    #
    shinyjs::reset('file_reference_gw_bi')
    shinyjs::reset('chr_select_bi')
    shinyjs::reset('start_interval_bi')
    shinyjs::reset('end_interval_bi')
    shinyjs::reset('run_bi')
    shinyjs::reset('bam_gw_bi')
    
    
    
    
  })
  
  observeEvent(input$tab == 'bisulfite', {
    
    #
    shinyjs::reset("file_reference")
    shinyjs::reset("files_problem")
    shinyjs::reset("fw_primer")
    shinyjs::reset("rv_primer")
    shinyjs::reset("primers")
    shinyjs::reset("num_tss")
    shinyjs::reset("tss")

    #
    shinyjs::reset('file_reference_gw')
    shinyjs::reset('chr_select')
    shinyjs::reset('start_interval')
    shinyjs::reset('end_interval')
    shinyjs::reset('run')
    shinyjs::reset('bam_gw')
    #
    shinyjs::reset('file_reference_gw_bi')
    shinyjs::reset('chr_select_bi')
    shinyjs::reset('start_interval_bi')
    shinyjs::reset('end_interval_bi')
    shinyjs::reset('run_bi')
    shinyjs::reset('bam_gw_bi')
    
    
  })
  
  observeEvent(input$tab == 'gw', {
    
    #
    shinyjs::reset("file_reference_bi")
    shinyjs::reset("files_problem_bi")
    #
    shinyjs::reset("file_reference")
    shinyjs::reset("files_problem")
    #
    shinyjs::reset('file_reference_gw_bi')
    shinyjs::reset('chr_select_bi')
    shinyjs::reset('start_interval_bi')
    shinyjs::reset('end_interval_bi')
    shinyjs::reset('run_bi')
    shinyjs::reset('bam_gw_bi')
    shinyjs::reset('file_vcf_bi')
    
    
  })
  
  observeEvent(input$tab == 'gw_bi', {
    
    #
    shinyjs::reset("file_reference_bi")
    shinyjs::reset("files_problem_bi")
    #
    shinyjs::reset("file_reference")
    shinyjs::reset("files_problem")
    #
    shinyjs::reset('file_reference_gw')
    shinyjs::reset('chr_select')
    shinyjs::reset('start_interval')
    shinyjs::reset('end_interval')
    shinyjs::reset('run')
    shinyjs::reset('bam_gw')
    shinyjs::reset('file_vcf')
    
  })
  
  
  #
  
  check_input_sanger <- function(input1, input2, input3, input4) {
    
    if (is.null(input1) & is.null(input2)) {
      "Please, upload a reference sequence"
    } else if (is.null(input3) & is.null(input4)) {
      
      "Please, upload problem sequences"
    } else {
      NULL
    }
  }

  
    # 1. Add reference sequence  
    
  f_reference <- reactive({
  
        if (input$tab == 'nome') {
          x <- input$file_reference
        } else if (input$tab == 'bisulfite') {
          x <- input$file_reference_bi
        }
    	
        data <- x$datapath
        file_ref <- read_seq_ref(data)
        CpGsites <- del_gcg_end(toString(file_ref), input$tab)
        GpCsites <- del_gcg_exo(toString(file_ref))
        result <- list("file_ref" = file_ref, "CpGsites" = CpGsites[[1]], "GpCsites" = GpCsites[[1]] )
        return(result)
  })
  
  
  # 2. Add problem sequences
  
  f_problem <- reactive({
    
    if (input$tab == 'nome') {
      x <- input$files_problem
    } else if (input$tab == 'bisulfite') {
      x <- input$files_problem_bi
    }
    
    validate(
      need(length(x) != 0, 'Please, use the left menu to upload your sequences')
    )
    
    seq_set <- DNAStringSet()
    store_names <- c()
    data_names_file <- c()
    if (length(x$datapath) != 1) {
     for (i in 1:nrow(x)) {
      data <- x$datapath[i]
      info_seq <- read_seq_pro(data)
      seq_set <- append(seq_set, DNAStringSet(info_seq$SEQ))
      data_names_file <- append(data_names_file , x$name[i])
      store_names <- append(store_names, substring(info_seq$NAME, 2))
     }
      
      seq_set@ranges@NAMES <- data_names_file
    } else {
      
      seq_set <- readDNAStringSet(x$datapath[1])
    }
    return(seq_set)
  })
  
  # 3. Quality Control
  
  qc_sequences <- reactive({

    qc_dataframe <- qc_seq(f_reference()$file_ref, f_problem(), convert = TRUE)
    return(qc_dataframe)
  })
  
  # 3. Trimming sequences with primers
  trim_sequences <- reactive({
    req(qc_sequences())
    req(input$primers == 1 | input$primers_bi == 1)
    
    if (input$tab == 'nome') {
      a <- input$fw_primer
      b <- input$rv_primer
      
    } else if (input$tab == 'bisulfite') {
      a <- input$fw_primer_bi
      b <- input$rv_primer_bi
      
    }
    
    validate(
      need(a != 'Forward primer bisulfite converted' &
             b != 'Reverse primer bisulfite converted',
           'Please, introduce a correct string of DNA')
    )
    
    validate(
      need(length(unique(str_split(toupper(a), '')[[1]] %in% c('A', 'T', 'C', 'G'))) == 1 &
             length(unique(str_split(toupper(b), '')[[1]] %in% c('A', 'T', 'C', 'G'))) == 1  ,
           'Do your primers have four of the possible nucleotides (A, T, C, G)? ')
    )
    
    validate(
      need(unique(str_split(toupper(a), '')[[1]] %in% c('A', 'T', 'C', 'G')) == TRUE &
             unique(str_split(toupper(b), '')[[1]] %in% c('A', 'T', 'C', 'G')) == TRUE  ,
           'Do your primers have four of the possible nucleotides (A, T, C, G)? ')
    )
    
    validate(
      need(nchar(a) > 1 & nchar(b) > 1, 
           'Do your primers have four of the possible nucleotides (A, T, C, G)? ')
    )
    

    result <- primer_seq(qc_sequences()[[2]], a, b)
    
    return(result)
  })
  
  # 3.5 Choose
  
  choose_option <- reactive({
    if (input$primers == 1 | input$primers_bi == 1) {
      result <- trim_sequences()
    }
    else {
      result <- qc_sequences()[[2]]
    }
    return(result)
  })
  
  # 4. Alignment of sequences
  
  align_sequences <- reactive({
    req(choose_option())
    seq_for_align <- qc_sequences()[[2]]
    if (input$primers == 1 | input$primers_bi == 1) {
      seq_for_align <- trim_sequences()
    }
    
    result_align <<- align_seq(f_reference()$file_ref, choose_option() , meth_end = TRUE, meth_exo = TRUE, discard_gcg = input$tab)
    result <- list("END" = result_align[[1]], "EXO" = result_align[[2]])
    return(result) 
  })
  
  # 5. Dataframe translation
  
df_translate <- reactive({
  req(align_sequences())
  
  if (input$tab == 'nome') {
    if (input$type_arrang == 1) {
      df_translating <-  translating_df(align_sequences()$EXO, size_wd = input$arrang_gpc , pos_nucleosome = TRUE)
    } else if (input$type_arrang == 2) {
      df_translating <-  translating_df(align_sequences()$END, size_wd = input$arrang_cpg , pos_nucleosome = TRUE)
    }
    
  } else if (input$tab == 'bisulfite') {
    df_translating <-  translating_df(align_sequences()$END, size_wd = input$wd_size_bi , pos_nucleosome = FALSE)
    
  }

  return(df_translating) 
  })
  
  # 6. Clustering
df_clustering <- reactive({
  req(df_translate())
  correct_order <- clustering_df(df_translate(), f_problem(), NULL)$order
  dendrogram_clustering <- clustering_df(df_translate(), f_problem(), TRUE)[[1]]
  final_df <- align_sequences()[[2]][correct_order,]
  df_end <- align_sequences()[[1]][correct_order,]
  df_translate_correct <- df_translate()[correct_order,]
  result <<- list("END" = df_end , "EXO" = final_df, "PLOT" = dendrogram_clustering,
                 "ORDER" = correct_order, "DF_TRANSLATE_CORRECT" = df_translate_correct)
  return(result)
})

  # 6.1 Clustering dataframe
  resume_df <- reactive({
    vector_names <<- f_problem()@ranges@NAMES
    vector_names_order <<- vector_names[as.numeric(df_clustering()$ORDER)]
    req(df_clustering())
    df <- data.frame("Initial order" = vector_names , 
                     "Clustering order" = vector_names_order)
    return(df)
  })

  # 7. Report

# 7.1 Lollipop graphic

tss_on <- reactive({
  if (input$tss == 2 & input$tss_bi == 2) {
   result <- NULL
} else {
  
  if (input$tab == 'nome') {
    result <- input$num_tss
  } else if (input$tab == 'bisulfite') {
    result <- input$num_tss_bi
    
  }
}
  return(result)
})

lines_on <- reactive({
  if (input$red_lines == 2) {
    result <- NULL
  }
  else {
    result <- TRUE
  }
  return(result)
})

df_lollipop <- reactive({

  
  validate(
    check_input_sanger(input$file_reference, input$file_reference_bi,
                       input$files_problem, input$files_problem_bi)
  )
  
  validate(
    check_lines_cpg(input$red_lines, input$type_arrang)
  )

  req(df_clustering())
  
   if (input$tss == 1 | input$tss_bi == 1) {
  
     validate(
       need(!is.na(input$num_tss) & !is.na(input$num_tss_bi),
            'Please, introduce a correct TSS')
     )
     
   }

  
  
  if (input$tab == 'nome') {
    result <- lollipop_svg_doble_lines4(df_clustering()[[1]], df_clustering()[[2]], df_clustering()[[5]],
                                        input$mean_nucleosome, tss_on(), lines_on(),
                                        input$size_draw, input$color_circles, input$color_lines)
  } else if (input$tab == 'bisulfite') {
    result <- lollipop_bi(df_clustering()[[1]], df_clustering()[[2]], df_clustering()[[5]],
                                        0.7, tss_on(), NULL,
                                        input$wd_size_bi, input$color_circles_bi, input$color_circles_bi)
  }

  return(result)
})


df_lollipop_svg <- reactive({
  req(df_clustering())
  png()
  if (input$tab == 'nome') {
    result <- lollipop_svg_doble_lines4(df_clustering()[[1]], df_clustering()[[2]], df_clustering()[[5]],
                                        input$mean_nucleosome, tss_on(), lines_on(),
                                        input$size_draw, input$color_circles, input$color_lines)
  } else if (input$tab == 'bisulfite') {
    result <- lollipop_bi(df_clustering()[[1]], df_clustering()[[2]], df_clustering()[[5]],
                          0.7, tss_on(), NULL,
                          input$wd_size_bi, input$color_circles_bi, input$color_circles_bi)
  }
  
  result_svg <- grid.export("lollipop_vectorial.svg")
  dev.off()
  return(result_svg)
})

# 7.2.1 Data Analyses (Endogenous methylation)

nucleoplot_plot1 <- reactive({
  req(align_sequences())
  result <- nucleoplot_plot(align_sequences()$END, tss = tss_on(), unmeth = FALSE, n_occupancy = FALSE)
  return(result)
})


# 7.2.2 Data Analyses (Exogenous methylation)

nucleoplot_plot2 <- reactive({
  req(align_sequences())
  result <- nucleoplot_plot(align_sequences()$EXO, tss = tss_on(), unmeth = FALSE, n_occupancy = TRUE)
  return(result)
})

nucleoplot_plot3 <- reactive({
  req(align_sequences())
  result <- nucleoplot_plot(align_sequences()$EXO, tss = tss_on(), unmeth = TRUE, n_occupancy = TRUE)
  return(result)
})

# GRAPHIC PLOT

output$plot_window_iowa <- renderPlot({

  validate(
    check_input_sanger(input$file_reference, input$file_reference_bi,
                       input$files_problem, input$files_problem_bi)

  )
  
  validate(
    check_lines_cpg(input$red_lines, input$type_arrang)
    
  )

   plot_window(df_clustering()$DF_TRANSLATE_CORRECT, input$arrang_gpc, f_problem()@ranges@NAMES, 
               df_clustering()$ORDER,
               tss_on(), "plot", align_sequences()$EXO)
  

})

plot_window_graphic <- reactive({
  req(df_clustering())
  plot_window(df_clustering()$DF_TRANSLATE_CORRECT, input$arrang_gpc, f_problem()@ranges@NAMES, 
              df_clustering()$ORDER,
              tss_on(), "plot", align_sequences()$EXO)

})

# GRAPHIC PLOT

output$plot_window_iowa_bi <- renderPlot({

    validate(
    check_input_sanger(input$file_reference, input$file_reference_bi,
                       input$files_problem, input$files_problem_bi)
  )

  plot_window_bi(df_clustering()$DF_TRANSLATE_CORRECT, input$wd_size_bi, 
                f_problem()@ranges@NAMES,
              tss_on(), "plot", align_sequences()$END, df_clustering()$ORDER)
  
})

plot_window_graphic_bi <- reactive({
  req(df_clustering())
  plot_window_bi(df_clustering()$DF_TRANSLATE_CORRECT, input$wd_size_bi, f_problem()@ranges@NAMES,
              tss_on(), "plot", align_sequences()$END, df_clustering()$ORDER) 
  
})



# GRAPHIC HEATMAP

output$heatmap_web <- renderPlot({
  

  validate(
    need(input$type_arrang == 1, 'For correct identification of nucleosome regions, please select clustering by GpC sites')
    )
  
  validate(
    check_input_sanger(input$file_reference, input$file_reference_bi,
                       input$files_problem, input$files_problem_bi)
  )
  
  validate(
    check_lines_cpg(input$red_lines, input$type_arrang)
  )
 
 
  plot_window(df_clustering()$DF_TRANSLATE_CORRECT,input$arrang_gpc,
              f_problem()@ranges@NAMES, df_clustering()$ORDER,
              tss_on(), "heatmap", align_sequences()$EXO)
  
})

heatmap_window_graphic <- reactive({
  req(df_clustering())
  
  plot_window(df_clustering()$DF_TRANSLATE_CORRECT,input$arrang_gpc,
              f_problem()@ranges@NAMES, df_clustering()$ORDER,
              tss_on(), "heatmap", align_sequences()$EXO)
  
})

# GRAPHIC HEATMAP

output$heatmap_web_bi <- renderPlot({
  validate(
    check_input_sanger(input$file_reference, input$file_reference_bi,
                       input$files_problem, input$files_problem_bi)
  )

  
  plot_window_bi(df_clustering()$DF_TRANSLATE_CORRECT, input$wd_size_bi, f_problem()@ranges@NAMES,
              tss_on(), "heatmap", align_sequences()$END, df_clustering()$ORDER)
  
})

heatmap_window_graphic_bi <- reactive({
  req(df_clustering())
  plot_window_bi(df_clustering()$DF_TRANSLATE_CORRECT, input$wd_size_bi, f_problem()@ranges@NAMES,
              tss_on(), "heatmap", align_sequences()$END, df_clustering()$ORDER)
  
})
# 2 - OUTPUT REACTIVE LOLLIPOP

output$instant_plot <- renderPlot(

  
  
  df_lollipop()
)

output$instant_plot_bi <- renderPlot(
  df_lollipop()
)

# CREATION OF DATAFRAME WITH PERCENTAGES
info_percentages <- reactive({
  req(align_sequences())
  result <- generation_table(align_sequences()$END, align_sequences()$EXO)
  return(result)
})

# MERGE OF TWO GRAPHICS

nucleoplot_plot4 <- reactive({
  req(align_sequences())
  req(info_percentages())
  result <- nomeplot_merge(info_percentages(), tss = tss_on())
  return(result)
})

# GRID OF GRAPHICS

nucleoplot_plot5 <- reactive({
  req(align_sequences())
  req(info_percentages())
  req(nucleoplot_plot4())
  result <- plot_grid(nucleoplot_plot1(), nucleoplot_plot3(), ncol = 2, 
                      nrow = 1)
  return(result)
})

nucleoplot_plot5_bi <- reactive({
  req(align_sequences())
  req(info_percentages())
  req(nucleoplot_plot4())
  result <- plot_grid(nucleoplot_plot1(), ncol = 1, 
                      nrow = 1)
  return(result)
})

# OUTPUT ZIP ARCHIVE OF RESULTS

output$sample_seq <- downloadHandler(
  
  filename = 'sample_nome_pcr.zip',
  
  content = function(result.zip) {
    
    fs <- paste0('./sample_seq/NOMe-PCR_module/',
                 list.files('sample_seq/NOMe-PCR_module'))
     zip(zipfile= result.zip, files=fs)
     contentType = "application/zip"
    
  })

output$sample_seq2 <- downloadHandler(
  
  filename = 'sample_nome_seq.zip',
  
  content = function(result.zip) {
    
    fs <- paste0('./sample_seq/NOMe-seq_module/',
                 list.files('sample_seq/NOMe-seq_module'))
    zip(zipfile= result.zip, files=fs)
    contentType = "application/zip"
    
  })

output$sample_seq3 <- downloadHandler(
  
  filename = 'sample_bisulfite_seq.zip',
  
  content = function(result.zip) {
    
    fs <- paste0('./sample_seq/Bisulfite_module/',list.files('sample_seq/Bisulfite_module'))
    zip(zipfile= result.zip, files=fs)
    contentType = "application/zip"
    
  })

output$downloadData_nome <- downloadHandler(

  filename = 'nome_analysis.zip',
  
  content = function(result.zip) {

    result1 <- write.csv(info_percentages()$END, "information_CpG.csv")
    result2 <- write.csv(info_percentages()$EXO, "information_GpC.csv")
    
    result3 <- render('nome_analysis.Rmd', html_pretty(), "nome_analysis.html")                
    fs <- c("information_CpG.csv", "information_GpC.csv",
            "nome_analysis.html")
    zip(zipfile= result.zip, files=fs)
    contentType = "application/zip"
  })

output$downloadData_nome_bi <- downloadHandler(
  
  filename = 'nome_analysis.zip',
  
  content = function(result.zip) {
    
    result1 <- write.csv(info_percentages()$END, "information_CpG.csv")
    
    result3 <- render('nome_analysis_bisulfite.Rmd', html_pretty(), "nome_analysis.html")                
    fs <- c("information_CpG.csv", 
            "nome_analysis.html")
    zip(zipfile= result.zip, files=fs)
    contentType = "application/zip"
  })


output$downloadData_graphic <- downloadHandler(
  filename = 'lollipop_vectorial.zip',
  
  content = function(result.zip) {
    

    result4 <- df_lollipop_svg()
    
    fs <- c("lollipop_vectorial.svg")
    zip(zipfile= result.zip, files=fs)
    contentType = "application/zip"
  })

output$downloadData_graphic_bi <- downloadHandler(
  filename = 'lollipop_vectorial.zip',
  
  content = function(result.zip) {
    
    
    result4 <- df_lollipop_svg()
    
    fs <- c("lollipop_vectorial.svg")
    zip(zipfile= result.zip, files=fs)
    contentType = "application/zip"
  })

output$downloadData_technical <- downloadHandler(
  filename = 'technical_analysis.zip',
  
  content = function(result.zip) {
    result1 <- write.csv(align_sequences()$END, "CpG_binary_dataframe.csv")
    result2 <- write.csv(align_sequences()$EXO, "GpC_binary_dataframe.csv")
    result3 <- render('technical_report.Rmd', html_pretty(), "technical_report.html")
    
    fs <- c("CpG_binary_dataframe.csv", "GpC_binary_dataframe.csv", "technical_report.html")
    zip(zipfile= result.zip, files=fs)
    contentType = "application/zip"
    
  })

output$downloadData_technical_bi <- downloadHandler(
  filename = 'technical_analysis.zip',
  
  content = function(result.zip) {
    
    result1 <- write.csv(align_sequences()$END, "CpG_binary_dataframe.csv")
    
    result3 <- render('technical_report_bisulfite.Rmd', html_pretty(), "technical_report.html")
    
    fs <- c("CpG_binary_dataframe.csv",  "technical_report.html")
    zip(zipfile= result.zip, files=fs)
    contentType = "application/zip"
  })

######################## GENOME-WIDE ###############################
######################## GENOME-WIDE ###############################
######################## GENOME-WIDE ###############################


# SELECT CHROMOSOME


observe({
  x <- list_chromosomes[[input$file_reference_gw]]
  updateSelectInput(session,"chr_select", choices = x)
})

observe({
  x <- list_chromosomes[[input$file_reference_gw_bi]]
  updateSelectInput(session,"chr_select_bi", choices = x)
})



# START INTERVAL 

start_interval <- reactive({
  
  req(input$run | input$run_bi)
  library(Rsamtools)

  if (input$tab == 'gw') {
    
    result <- isolate(as.numeric(gsub(",", "", input$start_interval)))
    
  } else if (input$tab == 'gw_bi') {
    
    result <- isolate(as.numeric(gsub(",", "", input$start_interval_bi)))
  }
  

  return(result)
})


# END INTERVAL

end_interval <- reactive({
  
  req(input$run | input$run_bi)
  
  if (input$tab == 'gw') {
    
    result <- isolate(as.numeric(gsub(",", "", input$end_interval)))
    
  } else if (input$tab == 'gw_bi') {
    
    result <- isolate(as.numeric(gsub(",", "", input$end_interval_bi)))
  }
  return(result)
})

# CHROMOSOME SELECTION

chr_select_iso <- reactive({
 
  req(input$run | input$run_bi)
  
  if (input$tab == 'gw') {
  
  if (input$chr_select != 'X' & input$chr_select != 'Y' ) {

    result <- isolate(as.numeric(input$chr_select))
  } else {
    result <- isolate(input$chr_select)
   }
  } else if (input$tab == 'gw_bi') {
    
    if (input$chr_select_bi != 'X' & input$chr_select_bi != 'Y' ) {
      
      result <- isolate(as.numeric(input$chr_select_bi))
    } else {
      result <- isolate(input$chr_select_bi)
    }
  }
  
  return(result)
})

# TSS GENOME WIDE

# tss_on_gw <- reactive({
#   if (input$tss_gw == 2) {
#     result <- 0
#   } else {
#     result <- input$num_tss_gw
#   }
#   return(result)
# })

reference_gw <- reactive({
    
    #req(input$tab == 'gw')
    input$run
    input$run_bi
    
    ## LIST CHROMOSOMES
    
    mm9 <- list('chr1' = 1, 'chr2' = 2,'chr3' = 3,'chr4' = 4,'chr5' = 5,'chr6' = 6,'chr7' = 7,'chr8' = 8,'chr9' = 9,'chr10' = 10,'chr11' = 11,'chr12' = 12,'chr13' = 13,
                'chr14' = 14,'chr15' = 15,'chr16' = 16,'chr17' = 17,'chr18' = 18,'chr19' = 19,'chrX' = 'X','chrY' = 'Y')
    
    if (input$tab == 'gw') {
      result <- isolate(getting_ref(input$file_reference_gw,chr_select_iso(), start_interval(), end_interval()))
      
    } else if (input$tab == 'gw_bi') {
      
      result <- isolate(getting_ref(input$file_reference_gw_bi,chr_select_iso(), start_interval(), 
                                    end_interval()))
    }
    
    # result <- isolate(getting_ref('mm9',chr_select_iso(), start_interval(), end_interval()))

    return(result)
    
  })


  sites_gw <- reactive({
    
    req(input$tab == 'gw' | input$tab == 'gw_bi')
    
    reference_input <- reference_gw()
    start_interval <- start_interval()
    test1 <<- reference_input
    test2 <<- start_interval
    test3 <<- chr_select_iso()
    
    if (input$snp_gw == 1) {
      
      if (input$tab == 'gw') {
        gpc_list <- unname(str_locate_all(reference_input, 'GC[A,T,C]'))[[1]][,2] - 1
        cpg_list <- unname(str_locate_all(reference_input, '[A,T,C]CG'))[[1]][,1] + 1
        
        initial_value_gpc <- as.numeric(gpc_list) + start_interval - 1
        initial_value_cpg <- as.numeric(cpg_list) + start_interval - 1
      } else {
        
        gpc_list <- unname(str_locate_all(reference_input, 'GC[A,T,C]'))[[1]][,2] -1
        cpg_list <- unname(str_locate_all(reference_input, 'CG'))[[1]][,1] +1
        
        initial_value_gpc <- as.numeric(gpc_list) + start_interval - 1
        initial_value_cpg <- as.numeric(cpg_list) + start_interval - 1
      }
    
    
    } else {
      
      library(VariantAnnotation)
      
      if (input$tab == 'gw') {
      pathway_vcf <- input$file_vcf
	    gpc_list <- unname(str_locate_all(reference_input, 'GC[A,T,C]'))[[1]][,2] -1
     	cpg_list <- unname(str_locate_all(reference_input, '[A,T,C]CG'))[[1]][,1] +1
    	initial_value_gpc <- as.numeric(gpc_list) + start_interval - 1
    	initial_value_cpg <- as.numeric(cpg_list) + start_interval - 1
      } else {

     	gpc_list <- unname(str_locate_all(reference_input, 'GC[A,T,C]'))[[1]][,2] -1
     	cpg_list <- unname(str_locate_all(reference_input, 'CG'))[[1]][,1] +1
    	
	    initial_value_gpc <- as.numeric(gpc_list) + start_interval - 1
    	initial_value_cpg <- as.numeric(cpg_list) + start_interval - 1
        pathway_vcf <- input$file_vcf_bi
      }
      pathway_vcf <- pathway_vcf$datapath
      
      validate(
        need(!is.null(pathway_vcf), 'Please, upload a VCF file')
      )

      
    
    gpc_list <- vcf_discard(pathway_vcf, chr_select_iso(), initial_value_gpc, gpc_list) 
    
    cpg_list <- vcf_discard(pathway_vcf, chr_select_iso(), initial_value_cpg, cpg_list) 
    
    } 
 
    result <- list("CPG" = cpg_list, "GPC" = gpc_list, "CPG2" = initial_value_cpg, "GPC2" = initial_value_gpc)
    
    return(result) 
})

 


 problem_gw_raw <- reactive({



   
   what <- c("rname", "strand", "pos", "qwidth", "seq")


   if (input$tab == 'gw') {
     
     bamPath <- input$bam_gw
     
   } else if (input$tab == 'gw_bi') {
     
     bamPath <- input$bam_gw_bi
   }
   
   bamPath <- bamPath$datapath

  bf <- BamFile(bamPath)
  
  chrom_identifiers <- seqnames(seqinfo(bf))[nchar(seqnames(seqinfo(bf))) < 4]
  
  if (grepl(chrom_identifiers[1], 'chr')) {
    
  which <- GRanges( seqnames = paste0('chr', chr_select_iso()), ranges =  IRanges(start_interval(), end_interval()))
  } else {
      which <- GRanges( seqnames = chr_select_iso(), ranges =  IRanges(start_interval(), end_interval()))
  }
   param <- ScanBamParam(which = which, what = what)
   
   
     problem_seq <- scanBam(file = bamPath, index = indexBam(bamPath), param = param)
     problem_seq <- collect_seq(problem_seq, bamPath,  start_interval(), end_interval(), chr_select_iso(), which)
     
   
   
   return(problem_seq)
  
})


gpc_ranges <- reactive({
  
  result <- data.frame('start' = problem_gw_raw()$pos, 'end'= (problem_gw_raw()$pos + problem_gw_raw()$qwidth), 
                           'width' = problem_gw_raw()$qwidth)

  vector_delete <<- which(between(result[,1], start_interval(), end_interval()) & between(result[,2], start_interval(), end_interval()))
                                                                                                                         
  result <- result[vector_delete,]
  result_list <- list("RANGES" = result, "DELETE" = vector_delete)
  return(result_list)

})
  
  
problem_gw <- reactive({
  
  result <- problem_gw_raw()
  # result <- result[gpc_ranges()[[2]]]
  result$pos <- result$pos[gpc_ranges()[[2]]]
  result$qwidth <- result$qwidth[gpc_ranges()[[2]]]
  result$seq <- result$seq[gpc_ranges()[[2]]]

  return(result)
  
})
  

binary_gpc <- reactive({
  
  # GET NEW GPC_RANGES
  
  gpc_ranges <- data.frame('start' = problem_gw()$pos, 'end'= (problem_gw()$pos + problem_gw()$qwidth), 
                       'width' = problem_gw()$qwidth)
  
  # COMMON THINGS (CpG and GpC)
  
  mat <- matrix(0, length(DNA_ALPHABET), length(DNA_ALPHABET))
  mat[1:4, 1:4] <- c(1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0,
                     0, 0, 1)
  rownames(mat) <- colnames(mat) <- DNA_ALPHABET[1:length(DNA_ALPHABET)]

  # CPG SITES

  cpg_gw <- matrix(NA, nrow = length(problem_gw()$seq), ncol = length(sites_gw()[[1]]))
  colnames(cpg_gw) <- sites_gw()[[1]]
  cpg_gw <- as.data.frame(cpg_gw)

  
  initial_value <- as.numeric(colnames(cpg_gw)) + start_interval() - 1

  for (i in 1:dim(cpg_gw)[1]) {

    cpg_gw[i,] <- ifelse(between(initial_value, 
                                 as.numeric(gpc_ranges[i,][1]), 
                                 as.numeric(gpc_ranges[i,][2])), 'NA', 2)
    
  }
  
  align_sequences <- pairwiseAlignment(problem_gw()$seq , reference_gw(), substitutionMatrix = mat,
                                       type = "local")
  

    for (i in 1:length(align_sequences)) {


    temp_problem <- Biostrings::pattern(align_sequences[i])
    temp_reference <- subject(align_sequences[i])
    temp_reference_start <- start(temp_reference)


    ref <- paste(temp_reference)
    str <- paste(temp_problem)
    gc_ref <- start(matchPattern('CG', ref))
    gc_str <- start(matchPattern('CG', str))

    if (length(gc_ref) == 0) next



    intersect_gc <- as.numeric(gc_ref %in% gc_str)

    if ('NA' %in% cpg_gw[i,]) {
      pos_change <- which(cpg_gw[i,] == 'NA')[1]

      cpg_gw[i,][pos_change:(pos_change+length(intersect_gc)-1)] <- intersect_gc
      cpg_gw[i,] <- ifelse(cpg_gw[i,] == 'NA', 2, cpg_gw[i,])
    }
    }
  


  # GPC SITES
  
 gpc_gw <- matrix(NA, nrow = length(problem_gw()$seq), ncol = length(sites_gw()[[2]]))
 colnames(gpc_gw) <- sites_gw()[[2]]
 gpc_gw <- as.data.frame(gpc_gw)

 initial_value <- as.numeric(colnames(gpc_gw)) + start_interval() - 1
 
 for (i in 1:dim(gpc_gw)[1]) {

    gpc_gw[i,] <- ifelse(between(initial_value, as.numeric(gpc_ranges[i,][1]), as.numeric(gpc_ranges[i,][2])), 'NA', 2)

 }  

align_sequences <- pairwiseAlignment(problem_gw()$seq , reference_gw(), substitutionMatrix = mat, 
                                     type = "local")

  for (i in 1:length(align_sequences)) {
  
  temp_problem <- Biostrings::pattern(align_sequences[i])
  temp_reference <- subject(align_sequences[i])
  temp_reference_start <- start(temp_reference)
  
  
  ref <- paste(temp_reference)
  str <- paste(temp_problem)
  gc_ref <- start(matchPattern('GC', ref))
  gc_str <- start(matchPattern('GC', str))
  
  intersect_gc <- as.numeric(gc_ref %in% gc_str)
  
  if (length(intersect_gc) == 0) {
    gpc_gw <- gpc_gw[-i,]
    next
  }
  
  if ('NA' %in% gpc_gw[i,]) {
  pos_change <- which(gpc_gw[i,] == 'NA')[1]
  
  gpc_gw[i,][pos_change:(pos_change+length(intersect_gc)-1)] <- intersect_gc
  gpc_gw[i,] <- ifelse(gpc_gw[i,] == 'NA', 2, gpc_gw[i,])
  }
  }

# Eliminate rows with 2 only CpG

vector_store <- c()
for (i in 1:nrow(cpg_gw)) {
  if(length(unique(as.vector(cpg_gw[i,] == 2))) == 1) {
    vector_store <- append(i, vector_store)
  }
}

if (length(vector_store) != 0) {

  cpg_gw <- cpg_gw[-vector_store,]
}

#  Eliminate rows with 2 only CpG
vector_store <- c()

for (i in 1:nrow(gpc_gw)) {
  if(length(unique(as.vector(gpc_gw[i,] == 2))) == 1) {
    vector_store <- append(i, vector_store)
  }
}

if (length(vector_store) != 0) {

  gpc_gw <- gpc_gw[-vector_store,]
}

validate(
  need(length(cpg_gw) != 0,
       'Please, select a region with more CpG sites')
)

validate(
  need(length(gpc_gw) != 0,
       'Please, select a region with more GpC sites')
)

# REMOVE ROWS WITH NAs

cpg_gw <- cpg_gw[rowSums(cpg_gw == 'NA') == 0,]
gpc_gw <- gpc_gw[rowSums(gpc_gw == 'NA') == 0,]
cpg_gw <- na.omit(cpg_gw)
gpc_gw <- na.omit(gpc_gw)


    
result = list(cpg_gw,gpc_gw) 
return(result)


})

df_translate_gw <- reactive({
  
  req(binary_gpc())
  
  if (input$tab == 'gw') {
    if (input$type_arrang_gw == 1) {
      df_translating <-  translating_df(binary_gpc()[[2]], size_wd = input$arrang_gpc_gw , pos_nucleosome = TRUE)
    } else if (input$type_arrang_gw == 2)   {
      df_translating <-  translating_df(binary_gpc()[[1]], size_wd = input$arrang_cpg_gw , pos_nucleosome = TRUE)
    }
  } else if (input$tab == 'gw_bi') {
      df_translating <-  translating_df(binary_gpc()[[1]], size_wd = input$arrang_cpg_gw_bi , pos_nucleosome = TRUE)
  }
  
  return(df_translating) 
})


df_clustering_gw <- reactive({
  req(df_translate_gw())
  correct_order <- clustering_df(df_translate_gw(), f_problem(), NULL)$order
  final_df <- binary_gpc()[[2]][correct_order,]
  df_end <- binary_gpc()[[1]][correct_order,]
  df_end <- na.omit(df_end)
  final_df <- na.omit(final_df)
  # df_translate_correct <- df_translate()[correct_order,]
  result <- list("END" = df_end , "EXO" = final_df)
  return(result)
})

# Lollipop graphics


df_lollipop_gw_gpc <- reactive({
  
  req(input$tab == 'gw')
  pushViewport(viewport(width = unit(1, 'npc'), height = unit(1, 'npc'), just = c(0.5,0.5), clip="off"))
  # result <- lollipop_gw(binary_gpc()[[2]], 2, 'blue', gpc = TRUE)
  result <- lollipop_gw(df_clustering_gw()[[2]], 2, 'blue', gpc = TRUE)
  # Por qué algunos funcionan sin min_interval de argumento y otras sí...
  
  return(result)
  
}) 

  df_lollipop_gw_cpg <- reactive({
    
    req(input$tab == 'gw')
    pushViewport(viewport(width = unit(1, 'npc'), height = unit(1, 'npc'), just = c(0.5,0.5), clip="off"))
    result <- lollipop_gw(df_clustering_gw()[[1]], 2, 'black', gpc = FALSE)

  
  return(result)

})
  
  df_lollipop_gw_cpg_bi <- reactive({
    
    req(input$tab == 'gw_bi')
    
    validate(
      check_interval(input$start_interval_bi, input$end_interval_bi)
      
    )
    
    validate(
      check_input(input$bam_gw, input$bam_gw_bi)
    )
    
    validate(
      
      need(input$run_bi, 'Click on the button: Run analysis')
    )
    
    validate(
      check_bam(problem_gw_raw()$seq)
    )
   
    
    pushViewport(viewport(width = unit(1, 'npc'), height = unit(1, 'npc'), just = c(0.5,0.5), clip="off"))
    result <- lollipop_gw(df_clustering_gw()[[1]], 2, 'black', gpc = FALSE)
    
    
    return(result)
    
  })
  
  
  check_interval <- function(input1, input2) {
    
     input1 <- as.numeric(str_remove_all(input1, ','))
     input2 <- as.numeric(str_remove_all(input2, ','))
   
    
    if (is.na(input1) | is.na(input2)) {
      "Please, select a correct interval on the left menu."
    } else if (input1 > input2) {
      "Please, select a correct interval on the left menu."
    } else if ((input2 - input1) > 2001 ) {
      "Please, select a correct interval with a maximum of 2000 nt."
    } else if (input1 == input2){
      "Please, select a correct interval on the left menu with different start and end."
    } else {
      NULL
    }
      
  }
  

  
  check_input <- function(input1, input2) {
    if (is.null(input1) & is.null(input2)) {
      "Please, upload a bam file and specify an interval with reads."
    } else {
      NULL
    }
  }
  
  check_lines_cpg <- function(input1, input2) {
    if (input1 == TRUE & input2 == 2) {
      "For correct identification of nucleosome regions, please select clustering by GpC sites."
    } else {
      NULL
    }
  }
  
  check_bam <- function(input1) {
    if (length(input1) == 0) {
      "Please, specify an interval with reads."
    } else {
      NULL
    }
  }

  
  
  
  
genomic_plot_prev <- reactive({
  
  # test4 <<- input$bam_gw
  # test4 <- test4$datapath
  # 
  # test5 <<- chr_select_iso()
  # test6 <<- sites_gw()[[3]]
  # test7 <<- sites_gw()[[4]]
  # test8 <<- start_interval()
  # test9 <<- end_interval()
  # test10 <<- input$file_reference_gw
  
  # plotviz('C:\\Users\\Requena\\Desktop\\sample_seq\\NOMe-seq_module\\NOMe-IMR90-Chr7-5314680-5822336.bam', test5, test6, test7,  
  #         test8, test9, FALSE, test10)
  
  

  if (input$tab == 'gw') {
    bamPath <- input$bam_gw
    bamPath <- bamPath$datapath
    result <- plotviz(bamPath, chr_select_iso(), sites_gw()[[3]], sites_gw()[[4]],  
                      start_interval(), end_interval(), FALSE, input$file_reference_gw)
    
  } else if (input$tab == 'gw_bi') {
    
    bamPath <- input$bam_gw_bi
    bamPath <- bamPath$datapath
    result <- plotviz(bamPath, chr_select_iso(), sites_gw()[[3]], sites_gw()[[4]],  
                      start_interval(), end_interval(), TRUE, input$file_reference_gw_bi)
  }

    return(result)

  })



  
# PLOTS OUTPUT
  
  
output$lollipop_gw_gpc <- renderPlot({
  
  req(input$tab == 'gw')
  
  validate(
    check_interval(input$start_interval, input$end_interval)
  )
  
  validate(
    check_input(input$bam_gw, input$bam_gw_bi)
  )
  
  validate(
    
    need(input$run, 'Click on the button: Start analysis')
  )
  
  validate(
    check_bam(problem_gw()$seq)
  )
  

  df_lollipop_gw_gpc()
  
})
  
output$lollipop_gw_cpg <- renderPlot({
  
  req(input$tab == 'gw')
  
  validate(
    check_interval(input$start_interval, input$end_interval)
    
  )
  
  validate(
    check_input(input$bam_gw, input$bam_gw_bi)
   
  )
  
  validate(
   
    need(input$run, 'Click on the button: Start analysis')
  )
  
  validate(
    
    check_bam(problem_gw()$seq)
  )
 
  df_lollipop_gw_cpg()
  
})  

output$lollipop_gw_cpg_bi <- renderPlot({
  df_lollipop_gw_cpg_bi()
  
})  

output$genomic_plot <- renderPlot({
  req(input$tab == 'gw')
  
  validate(
    check_interval(input$start_interval, input$end_interval)
    
  )
  
  validate(
    check_input(input$bam_gw, input$bam_gw_bi)
    
  )
  
  validate(
    
    need(input$run, 'Click on the button: Start analysis')
  )
  
  validate(
   
    check_bam(problem_gw()$seq)
  )
 
  genomic_plot_prev()

})

output$genomic_plot_bi <- renderPlot({
  req(input$tab == 'gw_bi')
  
  
  validate(
    check_interval(input$start_interval_bi, input$end_interval_bi)
    
  )
  
  validate(
    
    check_input(input$bam_gw, input$bam_gw_bi)
    
  )
  
  validate(
   
    need(input$run_bi, 'Click on the button: Start analysis')
  )
  
  validate(
    
    check_bam(problem_gw()$seq)
  )
  
  genomic_plot_prev()
  
})

output$plot_combined <- renderPlot({
  
  req(input$tab == 'gw')
  
  validate(
    check_interval(input$start_interval, input$end_interval)
    
  )
  
  validate(
    
    check_input(input$bam_gw, input$bam_gw_bi)
    
  )
  
  validate(
   
    need(input$run, 'Click on the button: Start analysis')
  )
  
  validate(
 
    check_bam(problem_gw()$seq)
  )
 
  
  together_gw(info_percentages_gw()$cpg, info_percentages_gw()$gpc)
  
})



output$plot_cpg_gw <- renderPlot({
  req(input$tab == 'gw')
  
  validate(
    check_interval(input$start_interval, input$end_interval)
   
  )
  
  validate(
    
    check_input(input$bam_gw, input$bam_gw_bi)
   
  )
  
  validate(
    
    need(input$run, 'Click on the button: Start analysis')
  )
  
  validate(
   
    check_bam(problem_gw()$seq)
  )
  
  
  graphics_gw(binary_gpc()[[1]],TRUE)

})

output$plot_cpg_gw_bi <- renderPlot({
  req(input$tab == 'gw_bi')
  
  validate(
    check_interval(input$start_interval_bi, input$end_interval_bi)
    
  )
  
  validate(
    
    check_input(input$bam_gw, input$bam_gw_bi)
  )
  
  validate(
    
    need(input$run_bi, 'Click on the button: Start analysis')
  )
  
  validate(
   
    check_bam(problem_gw()$seq)
  )
  graphics_gw(binary_gpc()[[1]],TRUE)
  
})
  

output$plot_gpc_gw <- renderPlot({
  req(input$tab == 'gw')
  

  validate(
    check_interval(input$start_interval, input$end_interval)
  )
  
  validate(
    
    check_input(input$bam_gw, input$bam_gw_bi)
  )
  
  validate(
    need(input$run, 'Click on the button: Start analysis')
  )
  
  validate(
    check_bam(problem_gw()$seq)
  )
  
  p <- graphics_gw(binary_gpc()[[2]], FALSE)
  
  p
})  

# DOWNLOAD SITES GENOME WIDE

info_percentages_gw <- reactive({
  
  result_cpg <- graphics_gw(binary_gpc()[[1]], cpg = TRUE, table = TRUE)
  result_gpc <- graphics_gw(binary_gpc()[[2]], cpg = FALSE, table = TRUE)
  result <- list('cpg' = result_cpg, 'gpc' = result_gpc)
  return(result)
})


df_lollipop_svg_gw_gpc <- reactive({
  
  png()
   result <- lollipop_gw(binary_gpc()[[2]], 2, 'blue', gpc = TRUE)

  # result_svg <- grid.export(name = "lollipop_vectorial_gpc.svg", addClasses = FALSE, progress = FALSE)
   result_svg <- grid.export(name = "lollipop_vectorial_gpc.svg")
  dev.off()
  return(result_svg)
})

df_lollipop_svg_gw_cpg <- reactive({
  
  
  result <- lollipop_gw(binary_gpc()[[1]], 2, 'black', gpc = FALSE)
  
  result_svg <- grid.export(name = "lollipop_vectorial_cpg.svg")
  dev.off()
  return(result_svg)
})

output$downloadData_graphic_gw <- downloadHandler(

  filename = 'lollipop_vectorial_gw.zip',

  content = function(result.zip) {

    result1 <- df_lollipop_svg_gw_gpc()
    result2 <- df_lollipop_svg_gw_cpg()


    fs <- c("lollipop_vectorial_gpc.svg", "lollipop_vectorial_cpg.svg")
    zip(zipfile= result.zip, files=fs)
    contentType = "application/zip"

  }) 

output$downloadData_graphic_gw_bi <- downloadHandler(
  
  filename = 'lollipop_vectorial_bi_gw.zip',
  
  content = function(result.zip) {
    
    result2 <- df_lollipop_svg_gw_cpg()
    
    
    fs <- c("lollipop_vectorial_cpg.svg")
    zip(zipfile= result.zip, files=fs)
    contentType = "application/zip"
    
  }) 



output$downloadData_gw <- downloadHandler(
  
  
  filename = 'nome_analysis.zip',
  
  content = function(result.zip) {
    
    result1 <- write.csv(info_percentages_gw()$cpg, "information_CpG.csv")
    result2 <- write.csv(info_percentages_gw()$gpc, "information_GpC.csv")
    
    result3 <- render('nome_analysis_gw.Rmd', html_pretty(), "nome_analysis_gw.html")                
    fs <- c("information_CpG.csv", "information_GpC.csv", "nome_analysis_gw.html")
    zip(zipfile= result.zip, files=fs)
    contentType = "application/zip"
  })

output$downloadData_gw_bi <- downloadHandler(
  
  
  filename = 'nome_analysis_bi.zip',
  
  content = function(result.zip) {
    
    result1 <- write.csv(info_percentages_gw()$cpg, "informationbi_CpG.csv")
    result3 <- render('bisulfite_gw.Rmd', html_pretty(), "bisulfite_gw.html")   
    
    fs <- c("informationbi_CpG.csv", "bisulfite_gw.html")
    zip(zipfile= result.zip, files=fs)
    contentType = "application/zip"
  })

})
