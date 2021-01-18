calculate_reduced_a_batch <- function(Zmat=NULL, varG=NULL, P=NULL, MMtsqrt=NULL, y=NULL, quiet=TRUE)
{

  ## internal function to AM
  if( !(nrow(P) ==  nrow(y))){
    message(" Error:  there is a problem with the  dimensions of  P, and/or the matrix y.")
    message("         They should have the same number of rows ")
    message(" The dimensions are: \n")
    message(" dim(P)      = ", dim(P), "\n")
    message(" dim(y)   = ", dim(y), "\n")
    return(NULL)

  }

 if(is.null(varG)){
   message(" VarG must be specified.")
   return(NULL)
   }

  if(is.null(P)){
   message(" P must be specified")
   return(NULL)
   }


  if(is.null(y)){
   message(" y must be specified")
   return(NULL)
   }
  
  if(is.null(Zmat)){
    a <- varG * MMtsqrt %*% P %*% y
  } else {
    a <- varG * MMtsqrt %*% t(Zmat) %*% P %*% y
  }

return(a)

}




