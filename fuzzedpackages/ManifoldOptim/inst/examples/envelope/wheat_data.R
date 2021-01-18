# This is based on the wheat protien example from: 
# 1.) Cook et. al. "Envelope Models for Parsimonious and Efficient 
# Multivariate Linear Regression" Statistica Sinica 2010.
# 2.) Cook et. al. "envlp: A MATLAB Toolbox for Computing Envelope Estimators
# in Multivariate Analysis
#
# website reference:
# https://www.jstatsoft.org/article/view/v062i08

require(MASS)
require(ldr)
library(Rcpp)
library(RcppArmadillo)
library(ManifoldOptim)
library(mvtnorm)

# start by loading in the data and building a data frame
wheat_data = read.table("wheatprotein.dat")

# fit an OLS linear model
m1 <- lm(cbind(V3, V4) ~ as.factor(V8), data=wheat_data)
summary(m1)

# look at estimates of beta
m1$coefficients[2,]

# calculate residual covariance
Sres <- cov(m1$residuals)

# calculate response marginal covariance
Sy <- cov(cbind(wheat_data$V3, wheat_data$V4))

# Notice the estimator for beta is:
# [7.522, -2.061]'
# with standard error of:
# [8.815, 9.681]'
# Conclude the standard error is too high to discern between the two groups.
# So, we need to use an envelope.

# assume envelope dimension is u = 1
n = dim(wheat_data)[1]; p = 1; r = 2; u = 1

# estimate model parameters and use grassman opt to get gamma_hat

# It looks like optimizer can return matrices that are not quite symmetric.
# So let's symmetrize before we evaluate.
tx <- function(x) {
	idx.gamma <- 1:r*u
	gammahat <- matrix(x[idx.gamma], r, u)
	list(gammahat = gammahat)
}

# objective function from equation 20
f <- function(x)
{
	# get Gammahat 
	par <- tx(x)
	Gammahat <- par$gammahat

	# calculate Gammahat_0 (SRM NOTE: only r-u evals & evects in Gammahat_0)
	tmp <- (diag(r)-Gammahat %*% t(Gammahat))
	evals <- eigen(tmp)

	# squareroot the evals, want Gammahat_0
	e_vals <- diag(sqrt(evals$values[1:(r-u)]), r-u)
	e_vecs <- matrix(evals$vectors[,1:(r-u)], nrow = r)
	Gammahat_0 = e_vecs %*% e_vals

	# The following should ALWAYS be an identity matrix
	#print(Gammahat %*% t(Gammahat) + Gammahat_0 %*% t(Gammahat_0))

	# calculate eq 20. of reference 1
	logD = log(det(t(Gammahat) %*% Sres %*% Gammahat)) + 
		log(det(t(Gammahat_0) %*% Sy %*% Gammahat_0))

	return(logD)

}

# Take the problem above (defined in R) and make a Problem object for it.
# The Problem can be passed to the optimizer in C++
mod <- Module("ManifoldOptim_module", PACKAGE = "ManifoldOptim")
prob <- new(mod$RProblem, f)

# initialize randomly
A = matrix(rnorm(r^2), nc=r); D = t(A)%*%A;
evals = eigen(D);
x0 = matrix(evals$vectors[,u],nr=r);

# Test the obj and grad fn
prob$objFun(x0)
#head(prob$gradFun(x0))

# optimize Gamma over Grassman
mani.params <- get.manifold.params(IsCheckParams = TRUE)
solver.params <- get.solver.params(DEBUG = 0, Tolerance = 1e-12,
	Max_Iteration = 1000, IsCheckParams = TRUE, IsCheckGradHess = FALSE)
mani.defn <- get.grassmann.defn(r,u)

res <- manifold.optim(prob, mani.defn, method = "LRBFGS", 
	mani.params = mani.params, solver.params = solver.params, 
	x0 = x0)

print(res)
print(tx(res$xopt))

# Calculate new beta estimate
result <- tx(res$xopt)

# Calculate P_sig_hat
P_sig_hat <- result$gammahat %*% t(result$gammahat) 

# Calculate beta_hat for envelope
beta_hat_env <- P_sig_hat %*% m1$coefficients[2,]
beta_hat_env

# Note that the output either matches page 12 of reference 2) exactly
#> beta_hat_env
#          [,1]
#[1,]  5.140503
#[2,] -4.678209

# or it's completely different depending on the initial random matrix
#> beta_hat_env
#         [,1]
#[1,] 2.404271
#[2,] 2.625680
