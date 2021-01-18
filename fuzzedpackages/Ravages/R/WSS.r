WSS <- function(x, genomic.region = x@snps$genomic.region) {

  if(!is.factor(genomic.region)) stop("'genomic.region' should be a factor")
  genomic.region <- droplevels(genomic.region)

  # MAF calculée comme dans le papier princeps
  Q <- (2*x@snps$N2+x@snps$N1 + 1) / ( 2*(x@snps$N0+x@snps$N1+x@snps$N2) +2 )
  W <- sqrt((nrow(x)-x@snps$NAs) * Q * (1-Q))
  weights.0 <- ifelse( Q > 0.5, 2/W, 0)
  weights.1 <- 1/W
  weights.2 <- ifelse( Q < 0.5, 2/W, 0)
    
  if(length(weights.0) != ncol(x) | length(weights.1) != ncol(x) | length(weights.2) != ncol(x) | length(genomic.region) != ncol(x)) {
    stop("x and weights dimensions mismatch")
  }
  genomic.region <- as.factor(genomic.region)
                            
  B <- .Call('oz_burden2', PACKAGE = "Ravages", x@bed, nlevels(genomic.region), genomic.region, weights.0, weights.1, weights.2)
  colnames(B) <- levels(genomic.region)
  rownames(B) <- x@ped$id
  return(B)
}


