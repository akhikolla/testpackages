library(glmnet)
library(SplitReg)
library(MASS)
context("Compare coefficients at output with glmnet for lambda_D = 0")
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

set.seed(1)
n <- 100
p <- 200
x_large <- matrix
rho = 0.8
sigma <- (1 - rho) * diag(x = 1, p, p) + rho
x_large <- mvrnorm(n, mu = rep(0, p), Sigma = sigma)
true_beta <- c(rep(1, 10), rep(0, p - 10))
y_large <- x_large %*% true_beta + rnorm(n)
x_large_std <- scale(x_large, scale = apply(x_large, 2, function(xj) { sqrt(mean((xj - mean(xj))^2))}))
y_large_cen <- y_large - mean(y_large)
y_large_std <- y_large_cen / sqrt(mean(y_large_cen**2))
# Groups to use for EN in all comparisons
num_groups <- 3
# alpha to use
alphas <- c(0, 1/2, 3/4, 1)

# for(alpha in alphas){
#   test_that(paste0("Equality for p<n, alpha = ", alpha), {
#     fit_glmnet <- glmnet(x_small_std, y_small_std, alpha = alpha)
#     
#     fit_split <- Ensemble_EN_Grid(x_small_std, y_small_std, which_lambda = 1, lambdas_grid = sort(fit_glmnet$lambda),
#                                 lambda_fixed = 0, alpha = alpha, num_groups = num_groups,
#                                 tolerance = 1e-10, max_iter = 1e10)
#     
#     coef_glmnet <- sapply(sort(fit_glmnet$lambda),
#                           function(k, fit_glmnet) { as.numeric(predict(fit_glmnet, type = 'coef', s = k))[-1]}, fit_glmnet)
#     coef_split <- sapply(1:length(fit_glmnet$lambda),
#                            function(k, fit_split) { apply(fit_split[,,k], 1, mean)}, fit_split)
#     
#     error <- (coef_glmnet - coef_split)^2
#     
#     
#     expect_lte(max(error), 1e-4)
#   })
#   
#   
#   test_that(paste0("Equality for p>n, alpha = ", alpha), {
# 
#     fit_glmnet <- glmnet(x_large_std, y_large_std, alpha = alpha)
#     
#     fit_phalanx <- Ensemble_EN_Grid(x_large_std, y_large_std, which_lambda = 1, lambdas_grid = sort(fit_glmnet$lambda),
#                                 lambda_fixed = 0, alpha = alpha, num_groups = num_groups,
#                                 tolerance = 1e-10, max_iter = 1e10)
#     
#     coef_glmnet <- sapply(sort(fit_glmnet$lambda),
#                           function(k, fit_glmnet) { as.numeric(predict(fit_glmnet, type = 'coef', s = k))[-1]}, fit_glmnet)
#     coef_phalanx <- sapply(1:length(fit_glmnet$lambda),
#                            function(k, fit_phalanx) { apply(fit_phalanx[,,k], 1, mean)}, fit_phalanx)
#     
#     error <- (coef_glmnet - coef_phalanx)^2
#     
#     
#     expect_lte(max(error), 1e-4)
#   })
# }