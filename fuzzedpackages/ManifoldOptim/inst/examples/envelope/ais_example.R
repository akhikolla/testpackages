# Date Created: 5/3/16
# Description:
# This is an implementation of an alternative algorithm based on product
# manifolds to algorithms 1 through 3 in:
# D. Cook, X. Zhang "Foundations for Envelope Models and Methods." 2014.
#
# Specifically this code uses the logistic regression example on the bottom 
# of page 23 on the Australia Institute of Sport (AIS) data set:
# http://www.statsci.org/data/oz/ais.html

require(MASS)
require(ldr)
library(Rcpp)
library(RcppArmadillo)
library(ManifoldOptim)
library(mvtnorm)

set.seed(1234)

# start by loading in the data and building a data frame
ais_data = read.table("ais.txt",header=TRUE)

# save off the number of data points
n <- nrow(ais_data)

# fit a logistic regression
m1 <- glm(as.factor(Sex) ~ Ht + Wt, data = ais_data, family="binomial")
summary(m1)

# look at estimates of beta
m1$coefficients

# setup dimensions
r <- 2
u <- 1

# marginal covariance of predictors
Sx <- cov(data.frame(Ht = ais_data$Ht, Wt = ais_data$Wt));

# backout gammahat, alphahat, and etahat
tx <- function(x) {
      
	idx.gamma <- 1:r*u
	idx.alpha <- 1:u + r*u
	idx.eta <- 1:u + u + r*u

	gammahat <- matrix(x[idx.gamma], r, u)
	alphahat <- matrix(x[idx.alpha], u, 1)
	etahat <- matrix(x[idx.eta], u, 1)
	list(gammahat = gammahat, alphahat = alphahat, etahat = etahat)

}

# calculate the objective function based on equation 3.7
f <- function(x) {
	par <- tx(x)
	gamma_hat <- par$gammahat
	alpha_hat <- par$alphahat
	eta_hat <- par$etahat

	# Calculate C_n
	C_n <- 0
	for( i in 1:n){
	
	     # back out Xi from Ht and Wt
	     xi <- matrix(as.numeric(ais_data[i,c(12,13)]), nrow=2)

	     # and back out Yi as a numeric value from the factor
	     yi <- (as.numeric(ais_data[i,1]) - 1)
	     #yi <- yi mod 2

	     # estimate beta
	     betahat <- gamma_hat %*% eta_hat;

	     # Calculate vi from page 14
	     vi <- alpha_hat + t(betahat) %*% xi;

	     # Calculate A(vi) from Table 1 on page 14
	     A_vi <- 1 + exp(vi);

	     # Calculate C(vi) from Table 1 on page 14
	     C_vi <- yi %*% vi - log(A_vi);

	     # if C_vi is infinite then skip
	     #if(is.infinite(C_vi)){
             #  C_n = C_n - 1e100
	     #  print("Infinite")
             #  next;
	     #}

	     # add this value to C_n
	     C_n = C_n + C_vi;
        }

	# Calculate M_n
	M_n <- -(n / 2) * (log(det(t(gamma_hat) %*% Sx %*% gamma_hat)) + log(det(t(gamma_hat) %*% solve(Sx) %*% gamma_hat)) + log(det(Sx)))

	#print(C_n)
	#print(M_n)

	# return the sum
	return(-(C_n + M_n))
}


# Take the problem above (defined in R) and make a Problem object for it.
# The Problem can be passed to the optimizer in C++
mod <- Module("ManifoldOptim_module", PACKAGE = "ManifoldOptim")
prob <- new(mod$RProblem, f)

# initialize randomly 
A = matrix(rnorm(r^2), nc=r); D = t(A)%*%A;
evals = eigen(D);
gamma0 <- matrix(evals$vectors[,1:u],nr=r);
alpha0 <- matrix(rnorm(u),nr=u);
eta0 <- matrix(rnorm(u),nr=u);
x0 <- c(gamma0, alpha0, eta0)

# Test the obj and grad fn
print(prob$objFun(x0))

mani.params <- get.manifold.params(IsCheckParams = TRUE)
solver.params <- get.solver.params(DEBUG = 2, Tolerance = 1e-300,
	Max_Iteration = 1000, IsCheckParams = TRUE, IsCheckGradHess = FALSE)
mani.defn <- get.product.defn(get.grassmann.defn(r,u), get.euclidean.defn(r,1), get.euclidean.defn(u,1))


res <- manifold.optim(prob, mani.defn, method = "LRBFGS", 
	mani.params = mani.params, solver.params = solver.params, x0 = x0)

#Max_Iteration = 1000
#
#method = "RTRSR1", 


print(res)

par <- tx(res$xopt) 
print(par)

print(par$gammahat %*% par$etahat)

# SRM NOTE: 
# using RTRSR1 algorithm, final obj value is 994.3577 and betahat_env is:
> print(par$gammahat %*% par$etahat)
           [,1]
[1,] 0.06250763
[2,] 0.09744533


# the final value for betahat_env in the text is (0.085, 0.086)'
# may be worth examining the standard error
