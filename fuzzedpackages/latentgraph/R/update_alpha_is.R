update_alpha_is <- function(lambda2, x_j, Xminusj_star, theta_hat, alpha_hat, delta_hat, Xj_star){
  # calulate eta, D and r
  eta <- Xminusj_star%*%theta_hat + Xj_star%*%alpha_hat + delta_hat
  D <- 1 / (1 + exp(-eta))
  r_alpha <- 4*x_j - Xminusj_star%*%theta_hat - delta_hat + eta - 4*D

  fit.alpha <- glmnet(Xj_star, r_alpha, standardize=FALSE, intercept=FALSE, lambda=lambda2, thresh=1e-20)
  alpha <- as.matrix(fit.alpha$beta)
  return(alpha)
}
