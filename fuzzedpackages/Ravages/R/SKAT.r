SKAT <- function(x, NullObject, genomic.region = x@snps$genomic.region,
                 weights = (1 - x@snps$maf)**24, maf.threshold = 0.5, 
                 get.moments = "size.based", estimation.pvalue = "kurtosis", 
                 params.sampling, cores = 10, debug = FALSE, verbose = TRUE){

  if(!is.factor(genomic.region)) stop("'genomic.region' should be a factor")
  genomic.region <- droplevels(genomic.region)

  #Test if phenotype is continuous or categorial
  if(NullObject$pheno.type == "categorial"){
    if(verbose) cat("Categorial phenotype \n")
    if(get.moments == "size.based"){
      get.moments <- NullObject$get.moments
    }else{
      if(!(get.moments %in% c("permutations", "bootstrap", "theoretical"))) stop("Wrong 'get.moments' specified")
    }

    if(verbose) cat(get.moments, "\n")
    
    if(get.moments == "bootstrap"){
      if(missing(params.sampling))
        params.sampling <- list(perm.target = 100, perm.max=5e4)
      res <- SKAT.bootstrap(x, NullObject, genomic.region = genomic.region,
                            weights = weights, maf.threshold = maf.threshold,
                            perm.target = params.sampling$perm.target, 
                            perm.max = params.sampling$perm.max, 
                            debug = debug, estimation.pvalue = estimation.pvalue)
    }
  
    if(get.moments == "permutations"){
      if(missing(params.sampling))
        params.sampling <- list(perm.target = 100, perm.max=5e4)
      res <- SKAT.permutations(x, NullObject, genomic.region = genomic.region,
                               weights = weights, maf.threshold = maf.threshold,
                               perm.target = params.sampling$perm.target, 
                               perm.max = params.sampling$perm.max, 
                               debug = debug, estimation.pvalue = estimation.pvalue)
    }

    if(get.moments == "theoretical"){
      if(missing(cores)) cores <- 10
      res <- SKAT.theoretical(x, NullObject, genomic.region = genomic.region, 
                              weights = weights, maf.threshold = maf.threshold,
                              estimation.pvalue = estimation.pvalue, cores = cores, debug = debug)
    }
  }
  if(NullObject$pheno.type == "continuous"){
    if(verbose) cat("Continuous phenotype \n")
    res <- SKAT.continuous(x, NullObject, genomic.region = genomic.region,
                           weights = weights, maf.threshold = maf.threshold,
                           estimation.pvalue = estimation.pvalue, debug = debug)
  }
  return(results = res)
}
