#====================================================================================================
# CALC_ARP - function for calculating ancestor relation posterior probabilities
#====================================================================================================
calc_arp <- function(data, max_parent_size){
  
  if (ncol(data) > 20){
    stop("Don't use on systems with more than 20 variables.")
  }

  # Zero-center data
  data <- apply(data, 2, function(x) x-mean(x))
  
  # Calculate parent scores and save as temporary files
  calc_parent_score_to_file(data, "fml", max_parent_size, file_out = "temp")
  
  # Calculate parent support using the APS solver
  aps_type <- "ar_modular"
  system(paste("aps-0.9.1/aps/aps", aps_type, "temp.fml.score temp.fml.arp", collapse = ""))
  
  # Read in calculated parent support from file
  ar_post <- as.matrix(read.delim("temp.fml.arp", header = FALSE, sep = " ", skip = ncol(data)+1))
  colnames(ar_post) <- NULL
  
  # Transform into normalized probabilities and transpose matrix
  ar_post <- t(exp(ar_post-ar_post[1,1]))
    
  # Delete temporary files
  file.remove("temp.fml.arp","temp.fml.score")
  
  return(ar_post)
}