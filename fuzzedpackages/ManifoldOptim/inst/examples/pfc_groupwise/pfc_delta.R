require(MASS)
require(ldr)
library(Rcpp)
library(RcppArmadillo)
library(ManifoldOptim)
library(mvtnorm)
source("MainFunctions_P.R")

set.seed(1234)

# ---- Start data generation ----

# start by constructing a structured covariance amoung predictors (i.e. groups exist)
n = 200; m = 4; p = 3*m
y <- array( rnorm(n=n, sd=3), c(n, 1));
fy <- scale(cbind(y, abs(y)), TRUE, FALSE)
A = matrix(rnorm(p^2), nc=p); D = t(A)%*%A
Gamma = eigen(D)$vectors[,1:2]

# inter-group covariances
rho <- c(.5, .5, .5)

# variances for all predictors
sigma2 <- 1

# group sizes
b <- c(4, 4, 4)

# build true Delta
tmp_idx = 1
Delta = matrix(0, nc=p, nr=p)

for(l in 1:length(b)){

ones = matrix(rep(1, b[l]), ncol=1)
Deltall = rho[l]*ones%*%t(ones) + (sigma2-rho[l])*diag(b[l])

#A = matrix(rnorm(b[l]^2,mean=0,sd=10),nc=b[l])
#Deltall = t(A)%*%A

Delta[tmp_idx:sum(b[1:l]), tmp_idx:sum(b[1:l])] <- Deltall
tmp_idx <- sum(b[1:l])+1

}




# Generate simulated data
Err <- t(mvrnorm(n=n, mu=c(rep(0, p)), Sigma=Delta));
X <- t(Gamma%*%t(fy) + Err); 



#--- End data generation ----



fy = bf(y, case="poly", degree=3) 

# calculate unstructured fit
fit0 = pfc(X, y, fy, numdir=2, structure="unstr")
fit0$Deltahat

# back out some additional values (SRM may be able to remove this...)
Xc <-  scale(X, TRUE, FALSE);
P_F <- fy%*%solve(t(fy)%*%fy)%*%t(fy);
Sigmahat <- cov(X)
Sfit <- t(Xc)%*%P_F%*%Xc/n; 
Sres <- Sigmahat - Sfit;



# It looks like optimizer can return matrices that are not quite symmetric.
# So let's symmetrize before we evaluate.

tx <- function(x) {
	S <- matrix(x, p, p)
	S[lower.tri(S)] <- t(S)[lower.tri(S)]
	return(S)
}


f <- function(x) {

	Deltahat <- tx(x)

	# calculate the likelihood associated with this candidate Delta
	# calculate eq 4. value for comparison
	#Deltahat <- fit0$Deltahat
	InvDeltahat <- solve(Deltahat)
	temp0 <- -(n*p/2)*log(2*pi);
	temp1 <- -(n/2)*log(det(Deltahat)); 
	temp2 <- -(n/2)*Trace(InvDeltahat%*%Sres);
	numdir = 2
	#if (numdir < p) temp3 <- -(n/2)*sum(eigen(InvDeltahat%*%Sres)$values[(numdir+1):p]);
	if (numdir < p)	temp3 <- - (n/2)*Trace(fit0$Gammahat)
	loglik <- Re(temp0 + temp1 + temp2 + temp3);
	return(-loglik)

}

# Take the problem above (defined in R) and make a Problem object for it.
# The Problem can be passed to the optimizer in C++
mod <- Module("ManifoldOptim_module", PACKAGE = "ManifoldOptim")
prob <- new(mod$RProblem, f)

mani.params <- get.manifold.params(IsCheckParams = TRUE)
solver.params <- get.solver.params(DEBUG = 0, Tolerance = 1e-12,
	Max_Iteration = 1000, IsCheckParams = TRUE, IsCheckGradHess = FALSE)
mani.defn <- get.spd.defn(p)


# -- optimization --
res <- manifold.optim(x0 = fit0$Deltahat, prob, method = "LRBFGS", mani.defn = mani.defn, mani.params = mani.params, solver.params = solver.params)


heatmap(Delta, Rowv=NA, Colv=NA)
heatmap(fit0$Deltahat, Rowv=NA, Colv=NA)
heatmap(tx(res$xopt), Rowv=NA, Colv=NA)


# --- when using sigma res and delta in term 3
#> f(cov(X))
#[1] 4176.678

#> f(Delta)
#[1] 4138.2

#> f(fit0$Deltahat)
#[1] 4096.354

# -- with cov(X) as initialization --
#> f(tx(res$xopt))
#[1] 4149.202

# -- with fit0$Deltahat as initialization --
#> f(tx(res$xopt))
#[1] 3761.271

# examine sum of eigen values in diff
Trace(fit0$Deltahat - Delta)
#[1] 0.3662055
Trace(tx(res$xopt) - Delta)
#[1] 10.37618


# --- when using gammahat in term 3

Trace(tx(res$xopt) - Delta)
# [1] 0.3731908

> f(fit0$Deltahat)
[1] 3081.955

> f(tx(res$xopt))
[1] 3081.332
