# checks the input object coefficients

stratEst.check.coefficients <- function( coefficients , covariates , num_strats , names_strategies ){

  if( is.matrix( coefficients ) ){
    rows_coefficients <- nrow( coefficients )
    cols_coefficients <- ncol( coefficients )

    if( rows_coefficients !=  length(covariates) ){
      stop("stratEst error: Input object 'coefficients' must have as many rows as there are covariates.");
    }
    if( cols_coefficients !=  num_strats ){
      stop("stratEst error: Input object 'coefficients' must have as many columns as there are strategies.");
    }
    if( any( is.na( coefficients ) ) ){
      stop("stratEst error: The input object 'coefficients' cannot contain NA values.");
    }
  }else{
    stop("stratEst error: Input object 'coefficients' has to be a matrix.");
  }

  coefficient_mat <- coefficients

  colnames_coefficients <- colnames(coefficient_mat)
  rownames_coefficients <- rownames(coefficient_mat)

  if( is.null(rownames_coefficients) == FALSE ){
    num.covariates <- length(covariates)
    names_covariates <- covariates
    new_coefficient_mat <- coefficient_mat
    for( v in 1:num.covariates){
      if( names_covariates[v] %in% rownames_coefficients  ){
        new_coefficient_mat[v,] <- coefficient_mat[names_covariates[v],]
      }else{
        stop(paste("stratEst error: There is a column named '", names_covariates[v] , "' in 'covariates' but no row with this name in 'coefficients'.",sep=""));
      }
    }
    coefficient_mat <- new_coefficient_mat
  }else{
    stop("stratEst error: The row names of the input object 'coefficients' must correspond to the column names of covariates.");
  }

  if( is.null(colnames_coefficients) == FALSE ){
    new_coefficient_mat <- coefficient_mat
    for( s in 1:num_strats){
      if( names_strategies[s] %in% colnames_coefficients  ){
        new_coefficient_mat[,s] <- coefficient_mat[,names_strategies[s]]
      }else{
        stop(paste("stratEst error: There is a strategy named '", names_strategies[s] , "' but no column with this name in 'coefficients'.",sep=""));
      }
    }
    coefficient_mat <- new_coefficient_mat
  }else{
    stop("stratEst error: The column names of the input object 'coefficients' must correspond to the names of the strategies.");
  }

  return( coefficient_mat )

}
