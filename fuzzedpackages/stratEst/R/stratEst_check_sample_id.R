# checks the input object sample.id

stratEst.check.sample.id <- function( data , sample.id ){

  if( is.null(sample.id) ){
    sample <- rep(1,nrow(data))
    num_samples <- 1
    sample_factor = FALSE
    sample_levels = FALSE
  }else{
    if( sample.id %in% colnames(data) ){
      sample <- data[,sample.id]
      if( "factor" %in% class(sample) == F ){
        warning(paste("stratEst error: The variable '",sample.id,"' ","has to be of class integer or factor.",sep=""))
        sample <- as.factor(sample)
      }
    }
    else{
      stop(paste("stratEst error: The data does not contain the variable '",sample.id,"' "," specified as sample id.",sep=""))
    }

    # sample
    sample_factor <- sample
    sample_levels <- as.character(unique(sample))
    # back to numeric
    sample <- as.numeric(sample)
    num_samples <- length(sample_levels)
    if( any( is.na( sample ) ) ){
      stop("stratEst error: The variable 'sample' in data cannot contain NA values.");
    }
  }

  stratEst.check.sample.id.return = list( "sample" = sample , "sample.factor" = sample_factor , "sample.levels" = sample_levels , "num.samples" = num_samples )

  return( stratEst.check.sample.id.return )

}
