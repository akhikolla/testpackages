burden.continuous <- function(x, NullObject, genomic.region = x@snps$genomic.region, burden, maf.threshold = 0.5){
  if(is.numeric(burden)) {
    if(!is.matrix(burden)){
      stop("Score is not a matrix")
    } else {
      if(is.null(colnames(burden))){ 
        colnames(burden) <- make.names(1:ncol(burden))
      }
    }
    score <- burden
  } else { 
    if(!is.factor(genomic.region)) stop("'genomic.region' should be a factor")
    genomic.region <- droplevels(genomic.region)
    if(missing(x)) stop("a bed.matrix 'x' is needed to compute the score")
    if(burden == "CAST"){
      score <- CAST(x, genomic.region, maf.threshold)
    } else if(burden == "WSS"){
      score <- WSS(x, genomic.region)
    } else {
      stop("'burden' should be \"CAST\", \"WSS\", or a matrix of pre-computed burdens");
    }
  }
  score <- as.data.frame(score)
  # to ensure syntactically correct formulas
  old.names <- colnames(score)
  names(score) <- make.names(names(score))

  # preparation data / formula
  data.reg <- cbind(NullObject$data, score) ; rownames(data.reg) <- NULL

  R <- sapply( names(score), function(reg) run.burden.continuous(pheno = NullObject$pheno, score = score, region = reg, covar.toinclude = NullObject$covar.toinclude, data = data.reg))

  R <- as.data.frame( t(R) );

  colnames(R) <- c("p.value", "is.err")

  rownames(R) <- old.names

  return(R)
}


run.burden.continuous <- function(pheno, score, region, covar.toinclude, data){
  # Formula for the current region
  if(is.null(covar.toinclude)) { 
    my.formula <- as.formula(paste("ind.pheno ~ ", region))
  } else {
    my.formula <- as.formula( paste("ind.pheno ~ ", region, " + ", covar.toinclude )) 
  }

 
  # Catch errors
  fit <- tryCatch(lm(my.formula, data = data), error = identity, warning = identity)

  if(is(fit, "error")) {
    pval <- NA ; 
    is.err <- 1 ; 
  } else {
    my.model <- summary(fit)
    pval <- my.model$coefficients[region, 4]
    is.err <- 0
  }

  return(c(pval, is.err))
}
  

  
