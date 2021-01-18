SKAT.NullObject <- function(group, data = NULL, formula){
  if(!is.factor(group))
    stop("'group' is not a factor")
    
  ref.level <- levels(group)[1]
  if(length(group)>=2000) get.moments <- "theoretical"

  if(is.null(data)) {
    a <- table(group)/length(group)
    Pi.data <- matrix( a, ncol = nlevels(group), nrow = length(group), byrow = TRUE)
    X <- matrix(1, nrow = length(group), ncol = 1) 
    if(length(group)<2000) get.moments <- "permutations"
  }else{
    if(!is.matrix(data)) stop("'data' should be a matrix")
    if(is.null(colnames(data))) colnames(data) <- sprintf("C%0*d", log10(ncol(data)) + 1, 1:ncol(data))
    if(is.null(formula)) formula <- as.formula(paste("~", paste(colnames(data), collapse = "+")))
    Pi.data <- Pi.matrix(group, data, formula, ref.level)
    X <- cbind(1, data[, all.vars(formula), drop=FALSE]) 
    if(length(group)<2000) get.moments <- "bootstrap"
  }
  #Compute matrix of var(y - ^pi)
  P1 <- P.mat2(Pi.data, X)  
  return(list(Pi.data = Pi.data, X = X, group = group, P1 = P1, get.moments = get.moments))
}


#For continuous phenotypes
SKAT.NullObject.continuous <- function(pheno, data = NULL, formula){
  if(is.null(data)){
    ymp <- matrix(lm(pheno ~ 1)$residuals, nrow = length(pheno), byrow = TRUE)
    X <- matrix(1, nrow = length(pheno), ncol = 1)
  }else{
    if(!is.matrix(data)) stop("'data' should be a matrix")
    if(is.null(colnames(data))) colnames(data) <- sprintf("C%0*d", log10(ncol(data)) + 1, 1:ncol(data))
    if(is.null(formula)) formula <- as.formula(paste("pheno ~", paste(colnames(data), collapse = "+")))
    #Residuals in regression directly gives y-pi_hat
    ymp <- matrix(lm(formula, data = as.data.frame(data))$residuals, ncol=1)
    #-1 to remove variable "pheno"
    X <- cbind(1, data[, all.vars(formula)[-1], drop=FALSE]) 
  }
  #Matrix P1
  V <- as.numeric(var(ymp)) * diag(length(pheno))
  P1 <- V - V %*% X %*% solve(t(X) %*% V %*% X) %*% t(X) %*% V

  return(list(ymp = ymp, X = X, pheno = pheno, P1 = P1))
}
