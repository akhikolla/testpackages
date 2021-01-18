## Generate a trivial dataset, X has mean 0 and norm 1, y has mean 0
set.seed(11)
n = 20
p = 5
x = matrix(rnorm(n*p), nrow=n, ncol=p)
x = scale(x, center = colMeans(x))
x = scale(x, scale = sqrt(colSums(x^2)))
beta = c(1, 1, 0, 0, 0)
y = x%*%beta + scale(rnorm(20, sd=0.01), center = TRUE, scale = FALSE)

## Perform 10-fold CV without replication
boss_cv_result = cv.boss(x, y)
## Get the coefficient vector of BOSS that gives minimum CV OSS score (S3 method for cv.boss)
beta_boss_cv = coef(boss_cv_result)
# the above is equivalent to
boss_result = boss_cv_result$boss
beta_boss_cv = boss_result$beta_boss[, boss_cv_result$i.min.boss, drop=FALSE]
## Get the fitted values of BOSS-CV (S3 method for cv.boss)
mu_boss_cv = predict(boss_cv_result, newx=x)
# the above is equivalent to
mu_boss_cv = cbind(1,x) %*% beta_boss_cv

## Get the coefficient vector of FS that gives minimum CV OSS score (S3 method for cv.boss)
beta_fs_cv = coef(boss_cv_result, method='fs')
## Get the fitted values of FS-CV (S3 method for cv.boss)
mu_fs_cv = predict(boss_cv_result, newx=x, method='fs')
