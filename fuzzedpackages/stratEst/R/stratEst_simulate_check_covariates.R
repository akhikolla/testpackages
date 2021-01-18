# checks the input object covariates

stratEst.simulate.check.covariates <- function( covariates ){

  if( "matrix" %in% class(covariates) == F ){
    stop(paste("stratEst error: The input object covariates must be a numeric matrix.",sep=""))
  }
  names_covariates <- colnames( covariates )
  if( is.null(names_covariates) ){
    stop(paste("stratEst error: The column names of the input object covariates have to indicate the names of the covariates.",sep=""))
  }
  covariate_mat <- covariates
  num_covariates <- ncol( covariate_mat )
  num_ids <- nrow(covariate_mat)

  stratEst.simulate.check.covariates.return <- list( "covariate.mat" = covariate_mat , "num.covariates" = num_covariates , "num.ids" = num_ids , "names.covariates" = names_covariates )

  return( stratEst.simulate.check.covariates.return )

}
