lambda.n = rev(exp(seq(1,45,1)/5 -7))
lambda.m = rev(exp(seq(1,45,1)/5 -7))
#lambda.n = rev(exp(seq(1,45,1)/4 -9))
#lambda.m = rev(exp(seq(1,45,1)/4 -9))
lambda.e = rev(exp(seq(1,45,1)/4 -9))
lambda.l = rev(exp(seq(1,45,1)/4 -9))

initiation <- function(x, y, alpha, family="gaussian"){
  lasso.cv <- glmnet::cv.glmnet(x,y, family=family, alpha=alpha, nfolds=5, intercept=FALSE)
  lambda <- lasso.cv$lambda.min
  lasso.fit <- glmnet::glmnet(x, y, family, alpha=alpha, nlambda=50, intercept=FALSE)
  coef0 <- as.vector(stats::predict(lasso.fit, s=lambda, type="coefficients", intercept=FALSE))[-1]
}

initiation_cox <- function(x, y0, d){
  y = cbind(time = y0, status = d)
  lasso.cv = glmnet::cv.glmnet(x, y, alpha=1, family="cox", nfolds=5, standardize=FALSE)
  alpha = 2*(lasso.cv$lambda.min)
  lasso.fit = glmnet::glmnet(x,y,family="cox", alpha=1, nlambda=100, standardize=FALSE)
  coef0 = as.numeric(stats::predict(lasso.fit, s=alpha, type="coefficients"))
}

TruePos <- function(b, b.true){
  index = which(b.true != 0)
  pos = which(b != 0)
  tp = length(intersect(index, pos))
  fp = length(pos) - tp
  list(tp=tp, fp=fp)
}

Adjacency = function(x, alpha=5)
{
  n = nrow(x)
  p = ncol(x)
  r0 = stats::cor(x)
  r = r0; r[which(r==1)] = 1 - 0.01
  z = 0.5*log((1+r[upper.tri(r)])/(1-r[upper.tri(r)]))
  c0 = mean(sqrt(n-3)*z) + 2*stats::sd(sqrt(n-3)*z)
  cutoff = (exp(2*c0/sqrt(n-3))-1)/(exp(2*c0/sqrt(n-3))+1)
  r = r0
  A = (r)^alpha*(abs(r)>cutoff)
  diag(A) = 0
  A
}

.onUnload <- function (libpath) {
  library.dynam.unload("regnet", libpath)
}
