## Generate a trivial dataset, X has mean 0 and norm 1, y has mean 0
set.seed(11)
n = 20
p = 5
x = matrix(rnorm(n*p), nrow=n, ncol=p)
x = scale(x, center = colMeans(x))
x = scale(x, scale = sqrt(colSums(x^2)))
beta = c(1, 1, 0, 0, 0)
y = x%*%beta + scale(rnorm(20, sd=0.01), center = TRUE, scale = FALSE)

## Fit the model
boss_result = boss(x, y)
## Print the values of AICc-hdf for all subsets given by BOSS
print(boss_result$IC_boss$aicc)
## calculate them manually using the calc.ic function
y_hat = cbind(rep(1,n),x)%*%boss_result$beta_boss
print(calc.ic(y_hat, y, df=boss_result$hdf_boss))
