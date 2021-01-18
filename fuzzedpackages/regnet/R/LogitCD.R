
LogitCD <- function(X, Y, penalty=c("network", "mcp", "lasso"), lamb.1=NULL, lamb.2=NULL, r=5, alpha=1,
                     init=NULL, alpha.i=1, standardize=TRUE)
{
  n = nrow(X); p = ncol(X);
  x = as.matrix(X); y = as.matrix(Y)
  b0 = rep(0, p+1)
  method = substr(penalty, 1, 1)
  #---------------------------------------------- Main Loop -----------------------------------------
  if(standardize) x = scale(x, scale = apply(x, 2, function(t) stats::sd(t)*sqrt((n-1)/n)))
  if(penalty == "network") a = Adjacency(x) else a = as.matrix(0)
  x = cbind(rep(1,n), x)
  init = match.arg(init, choices = c("zero","elnet"))
  if(init == "elnet") b0 = initiation(x, y, alpha.i, "binomial")

  b = RunLogit(x, y, lamb.1, lamb.2, b0, r, a, p, alpha, method)
  b = as.numeric(b)
  vname = colnames(x)
  if(!is.null(vname)){
    names(b) = c("Intercept", vname)
  }else{
    names(b) = c("Intercept", paste("v", seq = (1:p), sep=""))
  }

  return(drop(b))
}
