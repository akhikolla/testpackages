burden.weighted.matrix <- function(x, weights, genomic.region = x@snps$genomic.region){

  if(!is.factor(genomic.region)) stop("'genomic.region' should be a factor")
  genomic.region <- droplevels(genomic.region)

  if(length(weights)==1){
    weights <- rep(weights, length(genomic.region))
  }else{
    if(length(weights) != length(genomic.region)) stop("weights and genomic.region should have the same length")
  }

  weights.0 <- ifelse( x@snps$maf == (1 - x@p), 2*weights, 0)
  weights.1 <- weights
  weights.2 <- ifelse( x@snps$maf == x@p, 2*weights, 0)

  if(length(weights.0) != ncol(x) | length(weights.1) != ncol(x) | length(weights.2) != ncol(x) | length(genomic.region) != ncol(x)) {
    stop("x and weights dimensions mismatch")
  }
  genomic.region <- factor(genomic.region, levels = unique(genomic.region))

  B <- .Call('oz_burden2', PACKAGE = "Ravages", x@bed, nlevels(genomic.region), genomic.region, weights.0, weights.1, weights.2)
  colnames(B) <- levels(genomic.region)
  rownames(B) <- x@ped$id
  return(B)
}


