burden.NullObject <- function(group, ref.level, data, formula){
  if (!is.factor(group))  stop("'group' is not a factor")
  if (is.numeric(ref.level))  ref.level <- as.character(ref.level)
  if (!(ref.level %in% levels(group))) stop("'ref.level' is not a level of 'group'")
  if (is.null(data)){
        LogLik <- NULL
        covar.toinclude <- NULL
        data <- data.frame(ind.pheno = group) ; rownames(data) <- NULL
  }else{
    if(!is.matrix(data)) stop("'data' should be a matrix")
    if(is.null(colnames(data))) colnames(data) <- sprintf("C%0*d", log10(ncol(data))+1, 1:ncol(data))
    if (nrow(data) != length(group)) {stop("'data' has wrong dimensions")}
    if (is.null(formula)){
      covar.toinclude <- paste(colnames(data), collapse = "+")
      my.formula <- as.formula(paste("ind.pheno ~ 0 |", covar.toinclude))
    }else {
      z <- as.character(formula)
      if (z[1] != "~" | length(z) != 2)   stop("'formula' should be a formula of the form \"~ var1 + var2\"")
      covar.toinclude <- z[2]
      my.formula <- as.formula(paste("ind.pheno ~ 0 |", z[2]))
    }
    ##Creer data pour mlogit
    data <- as.data.frame(data) ; rownames(data) <- NULL
    data <- cbind(ind.pheno = group, data)
    data.reg <- dfidx(data, varying = NULL, shape = "wide", choice = "ind.pheno")
    ##Faire tourner le modele
    LogLik <- as.numeric(summary(mlogit(my.formula, data = data.reg, reflevel = ref.level))$logLik)
    
  }
  
  return(list(group = group, ref.level = ref.level, H0.LogLik = LogLik, covar.toinclude = covar.toinclude, data = data))
}


burden.NullObject.continuous <- function(pheno, data, formula){
  if (!is.numeric(pheno))  stop("'pheno' should be a numeric vector")
  if (is.null(data)){
        covar.toinclude <- NULL
        data <- data.frame(ind.pheno = pheno) ; rownames(data) <- NULL
  }else{
    if(!is.matrix(data)) stop("'data' should be a matrix")
    if(is.null(colnames(data))) colnames(data) <- sprintf("C%0*d", log10(ncol(data))+1, 1:ncol(data))
    if (nrow(data) != length(pheno)) {stop("'data' has wrong dimensions")}
    if (is.null(formula)){
      covar.toinclude <- paste(colnames(data), collapse = "+")
    }else {
      z <- as.character(formula)
      if (z[1] != "~" | length(z) != 2)   stop("'formula' should be a formula of the form \"~ var1 + var2\"")
      covar.toinclude <- z[2]
    }
    ##Creer data pour regression
    data <- as.data.frame(data) ; rownames(data) <- NULL
    data <- cbind(ind.pheno = pheno, data)
  }
  
  return(list(pheno = pheno, covar.toinclude = covar.toinclude, data = data))
}
