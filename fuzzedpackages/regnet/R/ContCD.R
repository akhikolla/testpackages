
ContCD <- function(X, Y, penalty=c("network", "mcp", "lasso"), lamb.1=NULL, lamb.2=NULL, clv=NULL, r=5, alpha=1,
                    init=NULL, alpha.i=1, standardize=TRUE)
{
  intercept = TRUE
  if(is.null(clv)){
    clv = intercept*1
  }else{
    clv = setdiff(union(intercept, (clv+intercept)), 0)
  }

  n = nrow(X); p.c = length(clv); p = ncol(X)-p.c+intercept;
  x = as.matrix(X); y = as.matrix(Y)
  b0 = rep(0, p+intercept)
  method = substr(penalty, 1, 1)
  #---------------------------------------------- Main Loop -----------------------------------------
  if(standardize) x = scale(x, scale = apply(x, 2, function(t) stats::sd(t)*sqrt((n-1)/n)))
  x = cbind(1, x)
  init = match.arg(init, choices = c("zero","elnet"))
  if(init == "elnet") b0 = initiation(x, y, alpha.i, "gaussian")

  x.c=x[, clv, drop = FALSE]; x.g = x[, -clv, drop = FALSE]
  if(penalty == "network") a = Adjacency(x.g) else a = as.matrix(0)

  b = RunCont(x.c, x.g, y, lamb.1, lamb.2, b0[clv], b0[-clv], r, a, p, p.c, method)
  # residual = y - cbind(x.c, x.g) %*% b
  b = as.numeric(b)
  vname = colnames(X)
  if(!is.null(vname)){
    names(b) = c("Intercept", vname[clv], vname[-clv])
  }else if(p.c==1){
    names(b) = c("Intercept", paste("g", seq = (1:p), sep=""))
  }else{
    names(b) = c("Intercept", paste("clv", seq = (1:(p.c-1)), sep=""), paste("g", seq = (1:p), sep=""))
  }

  # outlist = list(b=drop(b), residual=residual)
  # return(outlist)
  return(drop(b))
}

