update_alpha <- function(lambda2, x_j, Xminusj_star, theta_hat, delta_hat, Xj_star){
  r_alpha <- x_j - Xminusj_star%*%theta_hat - delta_hat
  fit.alpha <- glmnet(Xj_star, r_alpha, standardize=FALSE, intercept=FALSE, lambda=lambda2, thresh=1e-12)
  alpha <- as.matrix(fit.alpha$beta)
  return(alpha)
}
