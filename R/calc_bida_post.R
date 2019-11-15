#====================================================================================================
# Function for calculating the BIDA posteriors
#====================================================================================================
calc_bida_post <- function(data, ps){
  n <- nrow(data)
  d <- ncol(data)
  S <- t(data)%*%data
  nodes <- 1:d
  
  # Hyperparameters-------
  a0 <- 1
  b0 <- 1
  c_mu0 <- 0
  diag_V0 <- 1
  # ----------------------
  
  # Only include the most probable parent sets such that the cumulative support exceeds the given threshold
  th <- 0.999
  for (node in nodes){
    t <- sort(ps$parent_support[[node]], decreasing = TRUE, index.return = TRUE)		
    pos <- match(TRUE, cumsum(t$x) >= th)
    if (is.na(pos)){
      pos <- length(t$x)
    }
    ps$parent_support[[node]] <- t$x[1:pos]/sum(t$x[1:pos])
    ps$parent[[node]] <- ps$parent[[node]][t$ix[1:pos],]
  }
  
  # Calculate posteriors
  bida_post <- vector("list", length(d))	
  for (x in nodes){
    pars <- as.matrix(ps$parent[[x]])
    temp <- lapply(1:d, function(j) if(j != x){cbind(ps$parent_support[[x]],matrix(0,nrow(pars),3))})		
    for (i in 1:nrow(pars)){
      z <- pars[i,!is.na(pars[i,])]
      k <- length(z)
      mu0 <- matrix(c_mu0,k+1,1)
      A0 <- diag(k+1)/diag_V0
      xzTxz <- S[c(x,z),c(x,z)]
      A1 <- A0+xzTxz
      V1 <- solve(A1)
      a1 <- a0+n/2
      for (y in nodes[-x]){
        if (y %in% z) {
          temp[[y]][i,2:4] <- rep(0,1,3)
        } else {
          xzTy <- S[c(x,z),y]
          mu1 <- V1%*%(A0%*%mu0+xzTy)	
          b1 <- as.vector(b0+(t(mu0)%*%A0%*%mu0+S[y,y]-t(mu1)%*%A1%*%mu1)/2)
          v <- 2*a1
          V2 <- (b1/a1)*V1
          temp[[y]][i,2:4] <- c(mu1[1],V2[1,1],v)
        }				
      }			
    }
    for (y in nodes[-x]){
      ind <- which(temp[[y]][,3] == 0)
      if (length(ind) > 1){
        temp[[y]][ind[1],1] <- sum(temp[[y]][ind,1])
        temp[[y]] <- temp[[y]][-ind[-1],,drop = FALSE]		
      }	
    }
    bida_post[[x]] <- temp 		
  }		
  return(bida_post)		
}