
calculateMMt_sqrt_and_sqrtinv <- function(MMt=NULL, checkres=TRUE,
                                           quiet = TRUE )
{
  ## internal function to AM
  ## R function for calculating the square root of M * M^t
  ## and the inverse of the square root of MMt
  ## where M * M^t has already been created. 
  ## Using SVD for calculation
  ## Args
  ##  MMt   a matrix of M * M^t
  ##  checkres  when true, the accuracy of the inversion is checked. 
  
  ## testing that MMt is positive definite
  if(!is.positive.definite(MMt)){
    message(" Error: the matrix multiplication M %*% t(M) is not positive definite. \n")
    message("        This can occur if there are individuals with identical marker \n")
    message("        information. Please remove individuals with identical marker \n")
    message("        information, remembering also to remove their associated phenotype \n")
    message("        information as well. \n")
    message(" Internal function: calculateMMt_sqrt_and_sqrtinv has terminated with errors")
    return(NULL)
  } 
   res <- list()
      MMt.eigen <- eigen(MMt, symmetric=TRUE )
      sqrt_evals <- diag(sqrt(MMt.eigen$values))
      res[["sqrt"]] <- MMt.eigen$vectors %*% sqrt_evals %*% t(MMt.eigen$vectors)
      rm(MMt.eigen, sqrt_evals)
      gc()
      res[["invsqrt"]] <- chol2inv(chol(res[["sqrt"]]))



   if(checkres){
       a <- (res[["sqrt"]] %*% res[["invsqrt"]] )
       if(trunc(sum(diag(a))) != nrow(MMt))
       {
         message(" \n\n\nWARNING: these results may be unstable.\n")
         message(" The sum of the diagonal elements of the square root of M %*% t(M) and its inverse is ", sum(diag(a)), " where \n")
         message("  it should have been ", nrow(MMt), "\n")
         message("  This can occur if the genotype file contains near identical rows and/or columns.  Please check.\n\n")


       }
   }   ## end if(checkres)
   res <- list(sqrt_MMt=res[["sqrt"]], inverse_sqrt_MMt=res[["invsqrt"]] )



} ## end function



