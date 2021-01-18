update_deltaglm_is <- function(n,R,lambda3,p.fac, C_plus, x_j, Xminusj_star, theta_hat, Xj_star, alpha_hat, delta_hat){
  # calulate eta, D and r
  eta <- Xminusj_star%*%theta_hat + Xj_star%*%alpha_hat + delta_hat
  D <- 1 / (1 + exp(-eta))
  r_delta <- 4*x_j - Xminusj_star%*%theta_hat - Xj_star%*%alpha_hat + eta - 4*D

  fit.delta <- glmnet(x=C_plus, y=r_delta, standardize=FALSE, intercept=FALSE, lambda=lambda3/(n*R), penalty.factor=p.fac, thresh=1e-20)
  delta <- as.matrix(C_plus%*%fit.delta$beta)
  return(delta)
}
