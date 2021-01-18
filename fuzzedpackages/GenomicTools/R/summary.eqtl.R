summary.eqtl <- function(object, sig=0.01, ...){
  
  ws <- object$windowSize
  trans <- FALSE
  if(is.null(ws)){
    ws <- "trans-eQTL"
    trans <- TRUE
  }

  p.values <- sapply(object$eqtl,"[",3)
  
  if(trans){
    cat("trans-EQTL Summary\n")
  } else {
    cat("cis-EQTL Summary\n")
  }
  cat("---------------\n")
  cat("Type of test                                  :",object$method,"\n")

  testedGenes <- 0
  if(trans){
    testedGenes <- length(object$xAnnot)
  } else {
    testedGenes <- length(object$eqtl)
  }
  cat("Tested genes                                  :",testedGenes,"\n")

  cat("Total number of SNPs (in geno object)         :",dim(object$geno[[3]])[1],"\n")
  cat("Window size (in MB)                           :",ws,"\n")
  if(!trans) cat("Average (median) number of SNP in window      :",mean(sapply(p.values,length)),"(",median(sapply(p.values,length)),")\n")
  sigP <- sapply(p.values, function(x){sum(x<=sig, na.rm=TRUE)})
  if(!trans) cat("Average (median) number of sig. (p<",sig,") SNP in window :",mean(sigP),"(",median(sigP),")\n")
  invisible(object)
} 
