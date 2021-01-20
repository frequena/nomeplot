getting_ref <- function(genome_assembly, chrom, start, finish) {
  
  
  
  start <- as.character(start)
  finish <- as.character(finish)
  link <- paste0('http://genome.ucsc.edu/cgi-bin/das/',genome_assembly, '/dna?segment=', chrom, ':', start, ',', finish)
  
  r <- GET(link, content_type("text/plain"))
  stop_for_status(r)
  bin <- httr::content(r, "text")
  
  bin <- strsplit(bin, "\"")[[1]][17]
  bin <- gsub("[\r\n]", "", bin)
  bin <- gsub("[>]", "", bin)
  bin <- gsub("[<]", "", bin)
  bin <- gsub("[/]", "", bin)
  bin <- gsub("[A-Z]", "", bin)
  
  result <- DNAString(bin)
  return(result)
}




