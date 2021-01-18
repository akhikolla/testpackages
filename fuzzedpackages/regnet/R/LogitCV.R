
CV.Logit <- function(X, Y, penalty=c("network", "mcp", "lasso"), lamb.1=NULL, lamb.2=NULL, folds=5, r=5, alpha=1,
                    init=NULL, alpha.i=1, standardize=TRUE, verbo = FALSE)
{
  if(is.null(lamb.1)){
    lamb.1 = switch (penalty,
                     "network" = lambda.n,
                     "mcp" = lambda.m,
                     "lasso" = lambda.l)
  }
  if(is.null(lamb.2)) lamb.2 = c(0.1, 1, 10)
  init = match.arg(init, choices = c("zero","elnet"))

  n = nrow(X); p = ncol(X);
  X = as.matrix(X); Y = as.matrix(Y)

  b0 = rep(0, p+1)
  rs <- sample(c(1:n))
  CVM = matrix(0, length(lamb.1), length(lamb.2))
  method = substr(penalty, 1, 1)
  #---------------------------------------------- Main Loop -----------------------------------------
  for(f in 1:folds){
    if(verbo) cat("CrossValidation: ",f, "/", folds, "\n")
    index = c(1: ceiling(n/folds)) + (f-1)*ceiling(n/folds)
    test = rs[intersect(index, seq(1,n,1))]

    x = X[-test,]; y = Y[-test]
    x2 = X[test,]; y2 = Y[test]
    if(standardize){
      V1 = apply(x, 2, function(t) stats::sd(t)*sqrt((n-1)/n)); V1[V1==0]=1
      V2 = apply(x2, 2, function(t) stats::sd(t)*sqrt((n-1)/n)); V2[V2==0]=1
      x = scale(x, center = FALSE, scale = V1 )
      x2 = scale(x2, center = FALSE, scale = V2)
    }
    if(penalty == "network") a = Adjacency(x) else a = as.matrix(0)

    x = cbind(1, x); x2 = cbind(1, x2)
    if(init == "elnet") b0 = initiation(x, y, alpha.i, "binomial")

    CVM = CVM + LogitGrid(x, y, x2, y2, lamb.1, lamb.2, b0, r, a, p, alpha, method)

  }
  CVM = CVM/n
  mcvm = min(CVM)
  inds = which(CVM == mcvm, arr.ind=TRUE)
  lambda = lambda1 = lamb.1[inds[,1]]
  lambda2 = lamb.2[inds[,2]]
  if(length(lambda)>1) message("multiple optimal values(pairs) of lambda(s) are found.")
  rownames(CVM) = signif(lamb.1, digits = 3)
  if(penalty == "network"){
    lambda = cbind(lambda1, lambda2)
    colnames(CVM) = lamb.2
  }
  outlist = list(lambda=lambda, mcvm=mcvm, CVM=CVM, penalty=penalty)
  class(outlist) = "cv.logit"
  outlist
}
