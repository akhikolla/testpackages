update_theta <- function(lambda1, x_j, Xj_star, alpha_hat, delta_hat,Xminusj_star){
  r_theta <- x_j - Xj_star%*%alpha_hat - delta_hat
  fit.theta <- glmnet(Xminusj_star, r_theta, standardize=FALSE, intercept=FALSE, lambda=lambda1, thresh=1e-12)
  theta <- as.matrix(fit.theta$beta)
  return(theta)
}
