# This is a normal means comparison example based on the simulation 
# described on page 30 of:
# Cook et. al. "Envelope Models for Parsimonious and Efficient Multivariate 
# Linear Regression" 2009

require(MASS)
require(ldr)
library(Rcpp)
library(RcppArmadillo)
library(ManifoldOptim)
library(mvtnorm)

set.seed(1234)

n = 200; p = 1; r = 10; u = 1
beta = matrix(rep(sqrt(10),r),nr=r)

# Create a Gamma and Gamma_0
A = matrix(rnorm(r^2), nc=r); D = t(A)%*%A;
evals = eigen(D);
Gamma = matrix(evals$vectors[,u],nr=r);
Gamma_0 = evals$vectors[,(u+1):r];

# Create Omega and Omega_0
sigma = 3;
sigma_0 = 2;
Omega = sigma * diag(u);
Omega_0 = sigma_0 * diag(r-u);

# Create Sigma
Sigma = Gamma %*% Omega %*% t(Gamma) + Gamma_0 %*% Omega_0 %*% t(Gamma_0);

# Generate simulated data from population 1
y1 <- t(mvrnorm(n=n/2, mu=c(rep(0, r)), Sigma=Sigma));

# Generate simulated data from population 2
y2 <- t(mvrnorm(n=n/2, mu=beta, Sigma=Sigma));

# combine populations for total response
y <- t(cbind(y1,y2));

# combine predictors
x <- matrix(c(rep(0,n/2),rep(1,n/2)),nc=1)

# fit an OLS linear model
m1 <- lm(y ~ as.factor(x))
summary(m1)

# look at estimates of beta
m1$coefficients[2,]

# calculate residual covariance
Sres <- cov(m1$residuals)

# calculate response marginal covariance
Sy <- cov(y)

# It looks like optimizer can return matrices that are not quite symmetric.
# So let's symmetrize before we evaluate.
tx <- function(x) {
	idx.gamma <- 1:r*u

	gammahat <- matrix(x[idx.gamma], r, u)

	list(gammahat = gammahat)
}

# objective function from equation 20
f <- function(x) {

        # get Gammahat 
        par <- tx(x)
	Gammahat <- par$gammahat

	# calculate Gammahat_0 (SRM NOTE: only r-u evals & evects in Gammahat_0)
	tmp <- (diag(r)-Gammahat %*% t(Gammahat))
	evals <- eigen(tmp)

	# squareroot the evals, want Gammahat_0
	e_vals <- diag(sqrt(evals$values[1:(r-u)]), r-u) 
	e_vecs <- matrix(evals$vectors[,1:(r-u)],nr=r)
	Gammahat_0 = e_vecs %*% e_vals

	# The following should ALWAYS be an identity matrix
	#print(Gammahat %*% t(Gammahat) + Gammahat_0 %*% t(Gammahat_0))

	# calculate eq 20.
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

# determine the error for OLS
beta_hat_ols = matrix(m1$coefficients[2,],nr=r)
norm(beta_hat_ols - beta, "F")

#determine the error for Envelope
norm(beta_hat_env - beta, "F")
