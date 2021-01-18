library(glmnet)
library(SplitReg)
library(MASS)
context("Compare predict method with manual computation")
# Generate data sets, one with p<n, the other with p>n
set.seed(1)
n <- 100
p <- 20
x_small <- matrix
rho = 0.8
sigma <- (1 - rho) * diag(x = 1, p, p) + rho
x_small <- mvrnorm(n, mu = rep(0, p), Sigma = sigma)
true_beta <- c(rep(1, 10), rep(0, p - 10))
y_small <- x_small %*% true_beta + rnorm(n)
x_small_std <- scale(x_small, scale = apply(x_small, 2, function(xj) { sqrt(mean((xj - mean(xj))^2))}))
y_small_cen <- y_small - mean(y_small)
y_small_std <- y_small_cen / sqrt(mean(y_small_cen**2))
x_small_test <- x_small_std + rnorm(n)


set.seed(1)
n <- 50
p <- 75
x_large <- matrix
rho = 0.8
sigma <- (1 - rho) * diag(x = 1, p, p) + rho
x_large <- mvrnorm(n, mu = rep(0, p), Sigma = sigma)
true_beta <- c(rep(1, 10), rep(0, p - 10))
y_large <- x_large %*% true_beta + rnorm(n)
x_large_std <- scale(x_large, scale = apply(x_large, 2, function(xj) { sqrt(mean((xj - mean(xj))^2))}))
y_large_cen <- y_large - mean(y_large)
y_large_std <- y_large_cen / sqrt(mean(y_large_cen**2))
x_large_test <- x_large_std + rnorm(n)

# Groups to use for EN in all comparisons
groups <- c(2, 3, 5)


for(group in groups){
  test_that(paste0("Equality for p<n, groups = ", group), {
    object <- cv.SplitReg(x_small_std, y_small_std, num_models = group, num_folds = 5)
    preds <- predict(object, newx = x_small_test)
    coef <- object$betas[,,object$index_opt]
    preds_manual <- mean(y_small_std) - as.numeric(apply(apply(x_small_std, 2, mean) %*% coef, 1, mean)) + 
      as.numeric(apply(x_small_test %*% coef, 1, mean))
    preds_manual_cpp <- Prediction_Grid(x_small_test, x_small_std, y_small_std, object$betas)
    preds_manual_cpp <- apply(preds_manual_cpp[,,object$index_opt], 1, mean)
    error_R <- sqrt(sum((preds - preds_manual)^2))
    error_cpp <- sqrt(sum((preds - preds_manual_cpp)^2))
    expect_lte(max(error_R, error_cpp), 1e-10)
  })
  
  
  test_that(paste0("Equality for p>n, groups = ", group), {
    
    object <- cv.SplitReg(x_large_std, y_large_std, num_models = group, num_folds = 5)
    preds <- predict(object, newx = x_large_test)
    coef <- object$betas[,,object$index_opt]
    preds_manual <- mean(y_large_std) - as.numeric(apply(apply(x_large_std, 2, mean) %*% coef, 1, mean)) + 
      as.numeric(apply(x_large_test%*% coef, 1, mean))
    preds_manual_cpp <- Prediction_Grid(x_large_test, x_large_std, y_large_std, object$betas)
    preds_manual_cpp <- apply(preds_manual_cpp[,,object$index_opt], 1, mean)
    error_R <- sqrt(sum((preds - preds_manual)^2))
    error_cpp <- sqrt(sum((preds - preds_manual_cpp)^2))
    expect_lte(max(error_R, error_cpp), 1e-10)
  })
}