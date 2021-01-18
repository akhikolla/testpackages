update_theta_is <- function(lambda1, x_j, Xj_star, theta_hat, alpha_hat, delta_hat,Xminusj_star){
  # calulate eta, D and r
  eta <- Xminusj_star%*%theta_hat + Xj_star%*%alpha_hat + delta_hat
  D <- 1 / (1 + exp(-eta))
  r_theta <- 4*x_j - Xj_star%*%alpha_hat - delta_hat + eta - 4*D

  fit.theta <- glmnet(Xminusj_star, r_theta, standardize=FALSE, intercept=FALSE, lambda=lambda1, thresh=1e-20)
  theta <- as.matrix(fit.theta$beta)
  return(theta)
}
