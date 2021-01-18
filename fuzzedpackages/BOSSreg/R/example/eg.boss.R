## Generate a trivial dataset, X has mean 0 and norm 1, y has mean 0
set.seed(11)
n = 20
p = 5
x = matrix(rnorm(n*p), nrow=n, ncol=p)
x = scale(x, center = colMeans(x))
x = scale(x, scale = sqrt(colSums(x^2)))
beta = c(1, 1, 0, 0, 0)
y = x%*%beta + scale(rnorm(n, sd=0.01), center = TRUE, scale = FALSE)

## Fit the model
boss_result = boss(x, y)

## Get the coefficient vector selected by AICc-hdf (S3 method for boss)
beta_boss_aicc = coef(boss_result)
# the above is equivalent to the following
beta_boss_aicc = boss_result$beta_boss[, which.min(boss_result$IC_boss$aicc), drop=FALSE]
## Get the fitted values of BOSS-AICc-hdf (S3 method for boss)
mu_boss_aicc = predict(boss_result, newx=x)
# the above is equivalent to the following
mu_boss_aicc = cbind(1,x) %*% beta_boss_aicc

## Repeat the above process, but using Cp-hdf instead of AICc-hdf
## coefficient vector
beta_boss_cp = coef(boss_result, method.boss='cp')
beta_boss_cp = boss_result$beta_boss[, which.min(boss_result$IC_boss$cp), drop=FALSE]
## fitted values of BOSS-Cp-hdf
mu_boss_cp = predict(boss_result, newx=x, method.boss='cp')
mu_boss_cp = cbind(1,x) %*% beta_boss_cp
