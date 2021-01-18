update_deltaglm <- function(n,R,lambda3,p.fac, C_plus, x_j, Xminusj_star, theta_hat, Xj_star, alpha_hat){
  #n <- length(data)
  r_delta <- x_j - Xminusj_star%*%theta_hat - Xj_star%*%alpha_hat
  fit.delta <- glmnet(x=C_plus, y=r_delta, standardize=FALSE, intercept=FALSE, lambda=lambda3/(n*R), penalty.factor=p.fac, thresh=1e-12)
  delta <- as.matrix(C_plus%*%fit.delta$beta)
  return(delta)
}
