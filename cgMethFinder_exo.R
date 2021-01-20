cgMethFinder_exo <- function (ref, str) 
{
  ref <- paste(ref)
  str <- paste(str)
  cg_ref <- del_gcg_exo(ref)
  cg_str <- del_gcg_exo(str)
  indexcg_ref <- cg_ref[[1]][1:length(cg_ref[[1]])]
  indexcg_str <- cg_str[[1]][1:length(cg_str[[1]])]
  intersect_cg_cg <- as.numeric(indexcg_ref %in% indexcg_str)
  return(intersect_cg_cg)
}