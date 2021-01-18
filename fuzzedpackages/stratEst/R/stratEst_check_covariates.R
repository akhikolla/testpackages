# checks the input object covariates

stratEst.check.covariates <- function( data , covariates ){

  covariate_mat <- NULL
  for( i in 1:length(covariates) ){
    if( covariates[i] %in% colnames(data) ){
      covariate_vec <- data[,covariates[i]]
      if( "numeric" %in% class(covariate_vec) ){
        covariate_mat <- cbind(covariate_mat,covariate_vec)
      }else{
        stop(paste("stratEst error: The variable '",covariates[i],"' "," specified as covariate is not numeric.",sep=""))
      }
    }
    else{
      stop(paste("stratEst error: The data does not contain the variable '",covariates[i],"' "," specified as covariate.",sep=""))
    }
  }

  return( covariate_mat )

}
