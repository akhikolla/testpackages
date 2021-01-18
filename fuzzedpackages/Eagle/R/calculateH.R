calculateH <- function(MMt=NULL, varE=NULL, varG=NULL, Zmat=NULL)
{
  ## internal function to AM
  ## R function for calculating the H variance matrix 
  ## which is
  ##  H = \sigma^2_E I  + \sigma^2_G  MMt
  ## Args:
  ##     MMt  - matrix object for M %*% M^T
  ##     varE  -  numeric value for the residual variance
  ##     varG  -  numeric value for the polygenic variance (\sigma^2_g)
  ##
  ## H matrix is returned. 
  if(!is.numeric(varE)){
    message(" The varE (residual variance) must be numeric.")
    return(NULL)
    }

  if(varE < 0){
    message(" VarE cannot be negative.")
    return(NULL)
    }
  if(!is.numeric(varG)){
    message(" The varG (genotypic variance) must be numeric.")
    return(NULL)
    }
  if(varG < 0){
    message(" VarG cannot be negative.")
    return(NULL)
    }

  if(is.null(MMt)){
    message("MMt cannot be null.")
    return(NULL)
    }
  if(is.null(Zmat)){
  return( varE * diag(nrow(MMt)) + varG * MMt)
  } else {
  return( varE * diag(nrow(Zmat)) + varG * (Zmat %*% MMt %*% t(Zmat)) )

  }

}


