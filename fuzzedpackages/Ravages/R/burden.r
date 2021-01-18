burden <- function(x, NullObject, genomic.region = x@snps$genomic.region, burden, maf.threshold = 0.5, get.OR.value = FALSE, alpha = 0.05, cores = 10, verbose = TRUE){
  if(missing(x)) x <- NULL
  if(NullObject$pheno.type == "categorial"){
    if(verbose) cat("Categorial phenotype \n")
    res <- burden.mlogit(x = x, NullObject = NullObject, genomic.region = genomic.region, burden = burden, maf.threshold = maf.threshold, get.OR.value = get.OR.value, alpha = alpha, cores = cores)
  }
  if(NullObject$pheno.type == "continuous"){
    if(verbose) cat("Continuous phenotype \n")
    res <- burden.continuous(x = x, NullObject = NullObject, genomic.region = genomic.region, burden = burden, maf.threshold = maf.threshold)
  }
  return(res)
}
