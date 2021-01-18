# checks the input object coefficients

stratEst.simulate.check.coefficients <- function( coefficients , covariate_mat , num_strats , names_strategies ){

  if( is.matrix( coefficients ) ){
    rows_coefficients <- nrow( coefficients )
    cols_coefficients <- ncol( coefficients )

    if( rows_coefficients !=  ncol(covariate_mat) ){
      stop("stratEst error: Input object 'coefficients' must have as many rows as there are covariates.");
    }
    if( cols_coefficients !=  num_strats ){
      stop("stratEst error: Input object 'coefficients' must have as many columns as there are strategies.");
    }
    if( any( is.na( coefficients ) ) ){
      stop("stratEst error: The input object 'coefficients' cannot contain NA values.");
    }

    colnames_coefficients <- colnames(coefficients)
    rownames_coefficients <- rownames(coefficients)

    if( is.null(rownames_coefficients) == FALSE ){
      # if( "intercept" %in% rownames_coefficients  ){
      #   intercept <- coefficient_mat["intercept",]
      #   coefficient_mat <- rbind(intercept,coefficient_mat[rownames_coefficients!="intercept",])
      #   rownames_coefficients <- c("intercept",rownames_coefficients[rownames_coefficients!="intercept"])
      #   rownames(coefficient_mat) <- rownames_coefficients
      # }
      # else{
      #   stop("stratEst error: The input object 'coefficients' must contain a row named 'intercept'.");
      # }
      num.covariates <- ncol(covariate_mat)
      names_covariates <- colnames(covariate_mat)
      coefficient_mat <- coefficients
      for( v in 1:num.covariates){
        if( names_covariates[v] %in% rownames_coefficients  ){
          coefficient_mat[v,] <- coefficients[names_covariates[v],]
        }else{
          stop(paste("stratEst error: There is a column named '", names_covariates[v] , "' in 'covariates' but no row with this name in 'coefficients'.",sep=""));
        }
      }
    }else{
      stop("stratEst error: The row names of the input object 'coefficients' must correspond to the column names of covariates.");
    }

    if( is.null(colnames_coefficients) == FALSE ){
      for( s in 1:num_strats){
        if( names_strategies[s] %in% colnames_coefficients  ){
          coefficient_mat[,s] <- coefficients[,names_strategies[s]]
        }else{
          stop(paste("stratEst error: There is a strategy named '", names_strategies[s] , "' but no column with this name in 'coefficients'.",sep=""));
        }
      }
    }else{
      stop("stratEst error: The column names of the input object 'coefficients' must correspond to the names of the strategies.");
    }
  }else{
    stop("stratEst error: The input object 'coefficients' has to be a matrix.");
  }

  return( coefficient_mat )

}
