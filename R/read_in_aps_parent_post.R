#====================================================================================================
# Function for reading in APS computed parent posteriors
#====================================================================================================
read_in_aps_parent_posteriors <- function(file_path, max_par_size){
  temp <- read.table(file_path, header = FALSE, sep = " ", col.names = paste0("V",seq_len(max_par_size+2)), fill = TRUE)
  d <- temp[1,1]	
  pos <- 2
  parent <- vector("list", length = d)
  parent_support <- vector("list", length = d)
  for (i in 1:d){		
    node <- temp[pos,1]
    npar <- temp[pos,2]
    pos <- pos+1
    parent[[node]] <- temp[pos:(pos+npar-1),3:(max_par_size+2)]		
    s <- temp[pos:(pos+npar-1),1]
    s <- exp(s-max(s))
    parent_support[[node]] <- s/sum(s)		
    pos <- pos+npar		
  }
  ps <- list()
  ps$parent <- parent
  ps$parent_support <- parent_support
  return(ps)
}