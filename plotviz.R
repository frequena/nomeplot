

plotviz <- function(patway, chr_id, sites_genome_cpg, sites_genome_gpc,  start_id, end_id, only_cpg = FALSE,
                    genome_assembly) {
  
  import::from(Gviz, AnnotationTrack)
  import::from(Gviz, DataTrack)
  import::from(Gviz, plotTracks)
  import::from(Gviz, BiomartGeneRegionTrack)
  

  
  # test5 <<- chr_select_iso()
  # test6 <<- sites_gw()[[3]]
  # test7 <<- sites_gw()[[4]]
  # test8 <<- start_interval()
  # test9 <<- end_interval()
  # test10 <<- input$file_reference_gw
  
  # patway <- 'C:\\Users\\Requena\\Desktop\\sample_seq\\NOMe-seq_module\\Amh_2U.bam'
  # chr_id <- test5
  # sites_genome_cpg <- test6
  # sites_genome_gpc <- test7
  # start_id <- test8
  # end_id <- test9
  # genome_assembly <- test10
  # 
  # 
  # patway <- '/home/frequena/Desktop/archivos_BAM/in_chr3.bam'
  # chr_id <- manolo3
  # sites_genome_cpg <- manolo4
  # sites_genome_gpc <- manolo5
  # start_id <- manolo6
  # end_id <- manolo7
  # file.copy(patway, 'belero6.bam')
  # prueba <<- patway
  #file.copy(patway, 'temporal.bam')
  indexBam(patway)
  # patway <- paste0(patway, '.bam')
  
    # Configuration
  
  if(chr_id != 'X' & chr_id != 'Y') {
    chr_id <- as.numeric(chr_id)
  } else {
    chr_id <- paste0('chr', chr_id)
  }
  
  #patway <- paste0(patway, '.bam')
  #file.rename(patway, paste0(patway, '.bam'))
  # sites_genome_cpg <- as.character(sites_genome_cpg)
  # sites_genome_gpc <- as.character(sites_genome_gpc)
  
  plot_cpg <- GRanges(seqnames = chr_id, ranges = IRanges(sites_genome_cpg, sites_genome_cpg))
  plot_gpc <- GRanges(seqnames = chr_id, ranges = IRanges(sites_genome_gpc, sites_genome_gpc))
  # biomTrack <- BiomartGeneRegionTrack(genome= genome_assembly, chromosome=chr_id, start=start_id, 
  #                                     end=end_id, name="Ensembl", lwd = 6)
  
  # GGviz
  
  reads <-  Gviz::AnnotationTrack(range = patway, genome = genome_assembly, name = "Reads", chromosome = chr_id,  arrowHeadMaxWidth = 11)
  prueba1 <- Gviz::AnnotationTrack(range = plot_cpg, genome = genome_assembly, name = "CpG", chromosome = chr_id, fill = 'lightblue', fontsize = 10)
  prueba2 <- Gviz::AnnotationTrack(range = plot_gpc, genome = genome_assembly, name = "GpC", chromosome = chr_id, fill = 'red',  fontsize = 10)
  coverage <- Gviz::DataTrack(range = patway, genome = genome_assembly, type = "l", name = "Coverage", window = -1, chromosome = chr_id)
  
  
  if (only_cpg == FALSE) {
    list_plot <- list(prueba1, prueba2,reads, coverage)
  } else {
    
    list_plot <- list(prueba2,reads, coverage)
  }
  

  Gviz::plotTracks(list_plot, from = start_id, to = end_id, background.title = "steelblue", fontsize =20,
                   collapseTranscripts=T, shape="arrow", geneSymbols=T, arrowHeadWidth = 300)
  #file.remove('temporal.bam')
}


