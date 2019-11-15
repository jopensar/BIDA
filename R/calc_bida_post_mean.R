#====================================================================================================
# Function for calculating mean of BIDA posteriors
#====================================================================================================
calc_bida_post_mean <- function(bida_post){
  d <- length(bida_post)
  nodes <- 1:d	
  bida_mean <- matrix(0,d,d)	
  for (i in nodes){
    for (j in nodes[-i]){
      bida_mean[i,j] <- t(bida_post[[i]][[j]][,1])%*%bida_post[[i]][[j]][,2]			
    }
  }
  return(bida_mean)
}