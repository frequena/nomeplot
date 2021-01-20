# 
# 
#  patway <- "D:/nomeplot/bam_files/ACTB.bam"
# chr_id <- 7
#  sites_genome_cpg <- c(200, 300)
#  sites_genome_gpc <- c(200, 300)
#  start_id <- 5568253
#  end_id <- 5568626
#  genome_assembly <- 'hg19'


plotviz <- function(patway, chr_id, sites_genome_cpg, sites_genome_gpc,  start_id, end_id, only_cpg = FALSE,
                    genome_assembly) {
  
  import::from(Gviz, AnnotationTrack)
  import::from(Gviz, DataTrack)
  import::from(Gviz, plotTracks)
  import::from(Gviz, BiomartGeneRegionTrack)
  
  # 
  # manolo3 <<- input$chr_select
  # manolo4 <<- sites_gw()[[3]]
  # manolo5 <<-  sites_gw()[[4]]
  # manolo6 <<- start_interval()
  # manolo7 <<- end_interval()
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
  #patway <- paste0(patway, '.bam')
  
    # Configuration
  
  if(chr_id != 'X' & chr_id != 'Y') {
    chr_id <- as.numeric(chr_id)
  } else {
    chr_id <- paste0('chr', chr_id)
  }
  
  #patway <- paste0(patway, '.bam')
  #file.rename(patway, paste0(patway, '.bam'))
  sites_genome_cpg <- as.numeric(sites_genome_cpg)
  sites_genome_gpc <- as.numeric(sites_genome_gpc)
  
  plot_cpg <- GRanges(seqnames = chr_id, ranges = IRanges(sites_genome_cpg, sites_genome_cpg))
  plot_gpc <- GRanges(seqnames = chr_id, ranges = IRanges(sites_genome_gpc, sites_genome_gpc))
  biomTrack <- BiomartGeneRegionTrack(genome= genome_assembly, chromosome=chr_id, start=start_id, end=end_id, name="Ensembl", lwd = 6,
                                      filter = list(with_ox_refseq_mrna = TRUE))
  
  # GGviz
  
  reads <-  Gviz::AnnotationTrack(range = patway, genome = genome_assembly, name = "Reads", chromosome = chr_id,  arrowHeadMaxWidth = 11)
  prueba1 <- Gviz::AnnotationTrack(range = plot_cpg, genome = genome_assembly, name = "CpG", chromosome = chr_id, fill = 2, fontsize = 1,width = 10)
  prueba2 <- Gviz::AnnotationTrack(range = plot_gpc, genome = genome_assembly, name = "GpC", chromosome = chr_id, fill = 4,  fontsize = 1)
  coverage <- Gviz::DataTrack(range = patway, genome = genome_assembly, type = "l", name = "Coverage", window = -1, chromosome = chr_id)
  
  
  if (only_cpg == FALSE) {
    list_plot <- list(prueba1, prueba2, biomTrack,reads, coverage)
  } else {
    
    list_plot <- list(prueba2, biomTrack,reads, coverage)
  }
  
  Gviz::plotTracks(list_plot, from = start_id, to = end_id, background.title = "steelblue", fontsize =15,
                   collapseTranscripts=T, shape="arrow", geneSymbols=T, arrowHeadWidth = 300)
  #file.remove('temporal.bam')
}

# coverage <- Gviz::DataTrack(range = , genome = "mm9", type = "l", name = "Coverage", window = -1, chromosome = chr_id)
# Gviz::plotTracks(list(prueba1, prueba2, biomTrack, coverage), from = start_id, to = end_id, background.title = "steelblue", fontsize =15,
#                  collapseTranscripts=T, shape="arrow", geneSymbols=T, arrowHeadWidth = 300)
# BiomartGeneRegionTrack(genome="mm9", chromosome='X', start=83437700, end=83438570, name="Ensembl", lwd = 6)


