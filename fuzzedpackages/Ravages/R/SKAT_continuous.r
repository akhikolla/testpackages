SKAT.continuous <- function(x, NullObject, genomic.region = x@snps$genomic.region, 
                     weights = (1 - x@snps$maf)**24, maf.threshold = 0.5, 
                     estimation.pvalue = "kurtosis", debug = FALSE){

  if(!is.factor(genomic.region)) stop("'genomic.region' should be a factor")
  genomic.region <- droplevels(genomic.region)

  x@snps$weights <- weights
  x <- select.snps(x, x@snps$maf <= maf.threshold & x@snps$maf > 0)
  
  pheno <- matrix(NullObject$pheno, ncol = 1)
  
  #Matrix P1
  P1 <- NullObject$P1
  
  #(Y - pi_hat)
  ymp = NullObject$ymp
  
  #P-value for all regions
  res.allregions <- do.call(rbind, lapply(levels(genomic.region), function(reg) get.parameters.pvalue.continuous(x, region = reg, P1 = P1, ymp = ymp, estimation.pvalue = estimation.pvalue)))
  res.final <- as.data.frame(res.allregions)
  colnames(res.final) <- colnames(res.allregions)
  rownames(res.final) <- levels(genomic.region)
  if(debug)
    res.final
  else
    res.final[, c("stat", "p.value")]
} 
  

#P-value by genomic region
get.parameters.pvalue.continuous <- function(x, region, P1, ymp, estimation.pvalue){
  x.genomic.region <- select.snps(x, x@snps$genomic.region == region)
  GG <- gaston::as.matrix(x.genomic.region)
  #missing genotypes replaced by mean genotype
  GG <- apply(GG, 2, function(z) {z[which(is.na(z))] <- mean(z, na.rm = TRUE) ; return(z)})
  G <- GG %*% diag(x.genomic.region@snps$weights)
  
  #Stat de test
  Q <- as.vector(ymp) %*% (G %*% t(G)) %*% as.vector(ymp)

  #Moments
  M <- .Call("moments", PACKAGE = "Ravages", G, P1)

  #P-valeur
  pval <- p.valeur.moments.liu(Q = Q, mu = M["mu"], sigma = M["sigma"], skewness = M["skewness"], kurtosis = M["kurtosis"], estimation.pvalue = estimation.pvalue)

  return(c(stat = Q, p.value = pval, mean = as.numeric(M["mu"]), M["sigma"], M["skewness"], M["kurtosis"]+3))
}


  
