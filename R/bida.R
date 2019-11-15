#====================================================================================================
# BIDA main function
#====================================================================================================
bida <- function(data, max_parent_size){
  
  if (ncol(data) > 20){
    stop("Don't use on systems with more than 20 variables.")
  }

  # Zero-center data
  data <- apply(data, 2, function(x) x-mean(x))
  
  # Calculate parent scores and save as temporary files
  calc_parent_score_to_file(data, "fml", max_parent_size, fileOut = "temp")
  
  # Calculate parent support using the APS solver
  system(paste("aps-0.9.1/aps/aps", apsType,"temp.fml.score temp.fml.support", collapse = ""))
  
  # Read in calculated parent support from file
  ps <- read_in_aps_parent_post("temp.fml.support", max_parent_size)
  
  # Calculate BIDA posteriors
  bida_post <- calc_bida_post(data, ps)
  
  # Delete temporary files
  file.remove("temp.fml.support","temp.fml.score")
  
  return(bida_post)
}