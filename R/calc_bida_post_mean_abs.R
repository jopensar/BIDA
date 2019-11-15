#====================================================================================================
# Function for calculating (numerically) the mean absolute value of BIDA posteriors.
#====================================================================================================
calc_bida_post_mean_abs <- function(bida_post){
  d <- length(bida_post)
  nodes <- 1:d	
  bida_mean_abs <- matrix(0,d,d)
  df <- 0
  i <- 1
  while (df == 0){
    df <- max(sapply(bida_post[[i]][-i],function(x) max(x[,4])))
    i <- i+1
  }
  t <- seq(-5,5,by = 0.001)
  f <- dt(t,df)
  f <- f/sum(f)			
  for (i in nodes){
    for (j in nodes[-i]){
      temp <- bida_post[[i]][[j]]
      for (k in 1:nrow(temp)){
        x <- temp[k,2]+sqrt(temp[k,3])*t
        bida_mean_abs[i,j] <- bida_mean_abs[i,j]+temp[k,1]*sum(abs(x)*f)
      }
    }
  }
  return(bida_mean_abs)
}