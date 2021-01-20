


vcf_discard <- function(file, chromosome_number, first_sites, last_sites) {
  
  if (chromosome_number != 'X' | chromosome_number != 'Y' ) {
    
    chromosome_number <- as.numeric(chromosome_number)
  }
  prueba <- readVcf(file)
  prueba <- rowRanges(prueba)
  start(prueba)
  temp_chr <- as.numeric(as.vector(seqnames(prueba)))
  temp_logical <- unique(first_sites %in% start(prueba))
  
  if (length(temp_logical) == 2) {
    
    sites_delete_one <- first_sites[which(temp_chr == chromosome_number)]
    sites_delete_two <- first_sites[which(first_sites %in% start(prueba))]
    
    final_sites <- last_sites[!sites_delete_one %in% sites_delete_two]
    
  }
  
  return(final_sites)
  
}
