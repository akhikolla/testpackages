CAST <- function(x, genomic.region = x@snps$genomic.region, maf.threshold = 0.5) {

  if(!is.factor(genomic.region)) stop("'genomic.region' should be a factor")
  genomic.region <- droplevels(genomic.region)

  weights.0 <- ifelse( (x@snps$maf == (1-x@p) & x@snps$maf < maf.threshold), 1, 0)
  weights.1 <- ifelse( (x@snps$maf < maf.threshold), 1, 0)
  weights.2 <- ifelse( (x@snps$maf == x@p & x@snps$maf < maf.threshold), 1, 0)
      
  if(length(weights.0) != ncol(x) | length(weights.1) != ncol(x) | length(weights.2) != ncol(x) | length(genomic.region) != ncol(x)) {
    stop("x and weights dimensions mismatch")
  }
                
  B <- .Call('oz_burden2', PACKAGE = "Ravages", x@bed, nlevels(genomic.region), genomic.region, weights.0, weights.1, weights.2)
  colnames(B) <- levels(genomic.region)
  rownames(B) <- x@ped$id
  (B>0)+0
}

