# checks the input object shares

stratEst.check.shares <- function( shares , LCR , specific_shares , num_samples , num_strats , sample.id , sample_levels , select_strategies ){

  if( LCR ){
    shares = matrix( NA , num_strats , num_samples )
    warning("stratEst warning: Covariates specified. The input object 'shares' is ignored.");
  }else{
    if( specific_shares & num_samples > 1 ){
      if( inherits( shares , "list") == F ){
        stop("stratEst error: If the shares are sample specific, the input object 'shares' has to be a list of numeric vectors with as many elements as there are samples.");
      }
      shares_matrix = matrix( NA , num_strats , num_samples )
      sample_names_shares = names(shares)
      for( i in 1:num_samples ){
        expected_string = paste( sample.id , ".", sample_levels[i] , sep= "")
        if( expected_string %in% sample_names_shares ){
          shares_vec = shares[[expected_string]]
          if( length(shares_vec) != num_strats ){
            stop("stratEst error: The elements of the input object 'shares' have to be numeric vectors with as many elements as there are strategies.");
          }
          if( is.na(shares_vec[i]) == F & is.numeric(shares_vec[i]) == F  ){
            stop("stratEst error: The elements of the input object 'shares' have to be numeric vectors. NA values are allowed.");
          }
          shares_matrix[,i] = shares[[expected_string]]
        }else{
          stop(paste("stratEst error: There is no list element with name '", expected_string , "' in shares.", sep=""))
        }
      }
    }else{
      if( length(shares) != num_strats ){
        stop("stratEst error: The input object 'shares' has to be a numeric vector with as many elements as there are strategies.");
      }
      for( i in 1:num_strats ){
        if( is.na(shares[i]) == F & is.numeric(shares[i]) == F  ){
          stop("stratEst error: The input object 'shares' has to be a numeric vector. NA values are allowed.");
        }
      }
      shares_matrix = matrix( 0 , num_strats , num_samples )
      #for( i in 1:num_samples ){
        shares_matrix[,1] = shares
      #}
    }
    shares = shares_matrix
  }

  # check for fixed shares
  only_NA <- TRUE
  if( inherits( shares , "list") ){
    for( i in 1:length(shares)){
      if( any(is.na(shares[[i]]) == FALSE) ){
        only_NA = FALSE
      }
    }
  }else{
    if( specific_shares ){
      if( any(is.na(shares) == FALSE) ){
        only_NA = FALSE
      }
    }
    else{
      if( any(is.na(shares[,1]) == FALSE) ){
        only_NA = FALSE
      }
    }
  }

  if( select_strategies & only_NA == FALSE ){
    stop("stratEst error: If strategies are selected, shares cannot be fixed.");
  }

  return( shares )
}
