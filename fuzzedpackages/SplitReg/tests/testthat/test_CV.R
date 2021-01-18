library(glmnet)
library(SplitReg)
library(MASS)
context("Compare cross-validation over diversity and sparsity parameters")
# Generate data sets, one with p<n, the other with p>n
set.seed(1)
n <- 50
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
rho = 0.2
sigma <- (1 - rho) * diag(x = 1, p, p) + rho
x_large <- mvrnorm(n, mu = rep(0, p), Sigma = sigma)
true_beta <- c(rep(1, 10), rep(0, p - 10))
y_large <- x_large %*% true_beta + rnorm(n)
x_large_std <- scale(x_large, scale = apply(x_large, 2, function(xj) { sqrt(mean((xj - mean(xj))^2))}))
y_large_cen <- y_large - mean(y_large)
y_large_std <- y_large_cen / sqrt(mean(y_large_cen**2))
x_large_test <- x_large_std + rnorm(n)

num_groups <- 3
alpha <- 1
num_folds <- 5

ks <- c(1, 25, 35)
# test_that(paste0("Equality over sparsity for p<n"), {
#   set.seed(1)
#   lmax <- max(abs( t(x_small_std) %*% y_small_std / n))
#   lams <- exp(seq(log(1e-4 * lmax), log(lmax), length.out=50))
#   difs <- rep(NA, length(ks))
#   for (ind in 1:length(ks)){
#     k <- ks[ind]
#     cv_sparse <- CV_Ensemble_EN(x_small_std, y_small_std, 1, lams, 0,
#                                              alpha=alpha, num_groups=num_groups, num_folds=num_folds, tolerance=1e-11, max_iter=1e10,num_threads=1)
#     cv_div <- CV_Ensemble_EN(x_small_std, y_small_std,2, c(0,lams), lams[k],
#                                           alpha=alpha, num_groups=num_groups, num_folds=num_folds, tolerance=1e-11, max_iter=1e10,num_threads=1)
#     difs[ind] <- abs(1 - cv_div[1]/cv_sparse[k])
#   }
#   expect_lte(max(difs), 2e-2)
# })
# 
# test_that(paste0("Equality over diversity for p<n"), {
#   set.seed(1)
#   lmax <- max(abs( t(x_small_std) %*% y_small_std / n))
#   lams <- exp(seq(log(1e-4 * lmax), log(lmax), length.out=50))
#   difs <- rep(NA, length(ks))
#   for (ind in 1:length(ks)){
#     k <- ks[ind]
#     cv_sparse <- CV_Ensemble_EN(x_small_std, y_small_std, 1, lams, lams[1],
#                                 alpha=alpha, num_groups=num_groups, num_folds=num_folds, tolerance=1e-11, max_iter=1e10,num_threads=1)
#     cv_div <- CV_Ensemble_EN(x_small_std, y_small_std,2, c(0,lams), lams[k],
#                              alpha=alpha, num_groups=num_groups, num_folds=num_folds, tolerance=1e-11, max_iter=1e10,num_threads=1)
#     difs[ind] <- abs(1 - cv_sparse[k]/cv_div[2])
#   }
#   expect_lte(max(difs), 2e-2)
# })
# 
# test_that(paste0("Equality over sparsity for p>n"), {
#   set.seed(1)
#   lmax <- max(abs( t(x_large_std) %*% y_large_std / n))
#   lams <- exp(seq(log(1e-2 * lmax), log(lmax), length.out=50))
#   difs <- rep(NA, length(ks))
#   for (ind in 1:length(ks)){
#     k <- ks[ind]
#     cv_sparse <- CV_Ensemble_EN(x_large_std, y_large_std, 1, lams, 0,
#                                              alpha=alpha, num_groups=num_groups, num_folds=num_folds, tolerance=1e-11, max_iter=1e10,num_threads=1)
#     cv_div <- CV_Ensemble_EN(x_large_std, y_large_std,2, c(0,lams), lams[k],
#                                           alpha=alpha, num_groups=num_groups, num_folds=num_folds, tolerance=1e-11, max_iter=1e10,num_threads=1)
#     difs[ind] <- abs(1 - cv_div[1]/cv_sparse[k])
#   }
#   expect_lte(max(difs), 2e-2)
# })

# test_that(paste0("Equality over diversity for p>n"), {
#   set.seed(1)
#   lmax <- max(abs( t(x_large_std) %*% y_large_std / n))
#   lams <- exp(seq(log(1e-2 * lmax), log(lmax), length.out=50))
#   difs <- rep(NA, length(ks))
#   for (ind in 1:length(ks)){
#     k <- ks[ind]
#     cv_sparse <- CV_Ensemble_EN(x_large_std, y_large_std, 1, lams, lams[1],
#                                              alpha=alpha, num_groups=num_groups, num_folds=num_folds, tolerance=1e-11, max_iter=1e10,num_threads=1)
#     cv_div <- CV_Ensemble_EN(x_large_std, y_large_std,2, c(0,lams), lams[k],
#                                           alpha=alpha, num_groups=num_groups, num_folds=num_folds, tolerance=1e-11, max_iter=1e10,num_threads=1)
#     difs[ind] <- abs(1 - cv_sparse[k]/cv_div[2])
#   }
#   expect_lte(max(difs), 2e-2)
# })
