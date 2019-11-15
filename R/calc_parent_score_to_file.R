#====================================================================================================
# Calculates parent scores and writes to file in GOBNILP format. Assumes zero-centered data.
# Possible scores: 
# Fractional Marginal Likelihood ('fml') - hyperparameters: n0 = 1, a = d-1.
# Geiger & Heckerman score ('bge') - assumes a Wishart(I,d) prior.
#====================================================================================================
calc_parent_score_to_file <- function(data, score_type, max_parent_size, file_out){	

  # Zero-center data
  if (max(abs(apply(data, 2, mean))) > 1e-6){
    print("Assumes zero-centered data.")
    data <- apply(data, 2, function(x) x-mean(x))
  }
  
  n <- nrow(data)
  d <- ncol(data)
  S <- t(data)%*%data
  nodes <- 1:d	
  
  nps <- sum(sapply(0:max_parent_size, function(x) choose(d-1,x)))
  if (nps > 1e6){
    stop("The number of parent sets is > 10^6. Comment this out if you still want to proceed.")	
  }	
  
  fid <- file(paste(file_out, ".", score_type, ".score", sep = ""),"wt")
  writeLines(toString(d), con = fid, sep = "\n")
  if (score_type == "fml"){
    n0 <- 1
    a <- d-1
    const1 <- -log(pi)*(n-n0)*0.5
    for (node in nodes){
      writeLines(paste(node, nps), con = fid, sep = "\n")		
      for (k in 0:max_parent_size){
        dfam <- 1+k
        par_sets <- combn(nodes[-node],k)
        const2 <- lgamma(0.5*(a+n-d+dfam))-lgamma(0.5*(a+n0-d+dfam))+log(n0/n)*(a+n0-d+2*dfam-1)*0.5
        const <- const1+const2
        for (i in 1:ncol(par_sets)){				
          par <- par_sets[,i]
          fam <- c(node,par)
          scr <- const-(log(det(S[fam,fam, drop = FALSE]))-log(det(S[par,par, drop = FALSE])))*(n-n0)*0.5
          writeLines(paste(trimws(format(round(scr, 6), nsmall=6)),k,paste(par,collapse = " "), sep = " "), con = fid, sep = "\n")	
        }			
      }		
    }
    close(fid)
  } else if(score_type == "bge"){
    const1 <- -log(pi)*n*0.5
    for (node in nodes){
      writeLines(paste(node, nps), con = fid, sep = "\n")
      for (k in 0:max_parent_size){
        dfam <- 1+k
        par_sets <- combn(nodes[-node],k)
        const2 <- lgamma((dfam+n)*0.5)-lgamma(dfam*0.5)
        const <- const1+const2
        for (i in 1:ncol(par_sets)){
          par <- par_sets[,i]
          fam <- c(node,par)
          scr <- const-log(det(S[fam,fam, drop = FALSE]+diag(1,dfam)))*(dfam+n)*0.5+log(det(S[par,par, drop = FALSE]+diag(1,k)))*(k+n)*0.5
          writeLines(paste(trimws(format(round(scr, 6), nsmall=6)),k,paste(par,collapse = " "), sep = " "), con = fid, sep = "\n")                
        }
      }
    }
    close(fid)
  } else {
    file.remove(summary(fid)$description)
    close(fid)  
    stop("Unknown score type, possible options are 'fml' and 'bge'.")
  }
}