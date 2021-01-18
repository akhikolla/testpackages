# checks the input object shares

stratEst.simulate.check.shares <- function( shares , LCR , num_strats ){

  if( LCR ){
    shares = matrix( NA , num_strats , num_samples )
    warning("stratEst warning: Covariates specified. The input object 'shares' is ignored.");
  }else{
    if( "list" %in% class(shares) ){
      num_samples = length(shares)
      shares_matrix = matrix( NA , num_strats , num_samples )
      sample_names_shares = names(shares)
      null_names = is.null(sample_names_shares)
      for( i in 1:num_samples ){
        if( null_names ){
          sample_names_shares[i] = as.character(i)
          shares_vec = shares[[i]]
        }else{
          if( sample_names_shares[i] == "" ){
            sample_names_shares[i] = as.character(i)
            shares_vec = shares[[i]]
          }else{
            shares_vec = shares[[sample_names_shares[i]]]
          }
        }
        shares_matrix[,i] = shares_vec
      }
      names(shares) = sample_names_shares
    }else{
      num_samples = 1
      sample_names_shares = "1"
      shares_matrix = matrix( NA , num_strats , num_samples )
      shares_matrix[,1] = shares
      shares = shares_matrix
    }
  }

  stratEst.simulate.check.shares.return <- list("shares" = shares , "names_samples" = sample_names_shares , "num_samples" = num_samples )

  return( stratEst.simulate.check.shares.return )
}
