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

# desired reduction
numdirs = 2

# calculate unstructured fit
fit0 = pfc(X, y, fy, numdir=numdirs, structure="unstr")
fit0$Deltahat

# back out some additional values (SRM may be able to remove this...)
Xc <-  scale(X, TRUE, FALSE);
P_F <- fy%*%solve(t(fy)%*%fy)%*%t(fy);
Sigmahat <- cov(X)
Sfit <- t(Xc)%*%P_F%*%Xc/n; 
Sres <- Sigmahat - Sfit;



# It looks like optimizer can return matrices that are not quite symmetric.
# So let's symmetrize before we evaluate.
# SRM NOTE: also backing out two manifolds here (delta and gamma)
tx <- function(x) {
	idx.delta <- 1:p^2
	idx.gamma <- 1:p*numdirs + p^2


	deltahat <- matrix(x[idx.delta], p, p)
	deltahat = (deltahat + t(deltahat)) / 2
	gammahat <- matrix(x[idx.gamma], p, numdirs)

	list(deltahat = deltahat, gammahat = gammahat)
}

f <- function(x) {
        par <- tx(x)

	Deltahat <- par$deltahat
	Gammahat <- par$gammahat

	# calculate the likelihood associated with this candidate Delta
	# and candidate gamma
	# calculate eq 3. value for comparison

	InvDeltahat <- solve(Deltahat)
	InvDeltahatSqrt <- chol(InvDeltahat)
	temp0 <- -(n*p/2)*log(2*pi);
	temp1 <- -(n/2)*log(det(Deltahat)); 

	#print(t(Gammahat) %*% InvDeltahatSqrt %*% (InvDeltahatSqrt %*% Sfit %*% InvDeltahatSqrt) %*% InvDeltahatSqrt %*% Gammahat)
	#temp2 <- -(n/2)*Trace(InvDeltahatSqrt %*% Sigmahat %*% InvDeltahatSqrt -  InvDeltahatSqrt %*% Gammahat %*% t(Gammahat) %*% InvDeltahatSqrt %*% (InvDeltahatSqrt %*% Sfit %*% InvDeltahatSqrt))

	#temp2 <- -(n/2)*Trace(InvDeltahat %*% Sres) - (n/2)*Trace(Gammahat)

	# Kofi's equation
	temp2 <- -(n/2)*Trace(Sigmahat %*% InvDeltahat - t(Gammahat) %*% InvDeltahat %*% Sfit %*% Gammahat %*% solve(t(Gammahat) %*% InvDeltahat %*% Gammahat))

	loglik <- Re(temp0 + temp1 + temp2);
	return(-loglik)

}



# Take the problem above (defined in R) and make a Problem object for it.
# The Problem can be passed to the optimizer in C++
mod <- Module("ManifoldOptim_module", PACKAGE = "ManifoldOptim")
prob <- new(mod$RProblem, f)

# initialize with PFC output
x0 = c(as.numeric(fit0$Deltahat), as.numeric(fit0$Gammahat))

# Test the obj and grad fn
prob$objFun(x0)
#head(prob$gradFun(x0))

mani.params <- get.manifold.params(IsCheckParams = TRUE)
solver.params <- get.solver.params(DEBUG = 0, Tolerance = 1e-12,
	Max_Iteration = 1000, IsCheckParams = TRUE, IsCheckGradHess = FALSE)
mani.defn <- get.product.defn(get.spd.defn(p), get.grassmann.defn(p,2))

res <- manifold.optim(prob, mani.defn, method = "LRBFGS", 
	mani.params = mani.params, solver.params = solver.params, 
	x0 = x0)

print(res)
print(tx(res$xopt))


# examine likelihood
prob$objFun(res$xopt)

#[1] 2949.486



# examine heatmap
heatmap(tx(res$xopt)$deltahat, Rowv=NA, Colv=NA)

# examine sum of eigen values in diff
Trace(fit0$Deltahat - Delta)
#[1] 0.3662055
Trace(tx(res$xopt)$deltahat - Delta)
#[1] 0.4785665

# examine the sum of singular values in diff
sum(svd(Gamma - fit0$Gammahat)$d)
#[1] 3.994087
sum(svd(Gamma - tx(res$xopt)$gammahat)$d)
#[1] 3.07502



# calculate eq 4. value for comparison
Deltahat = tx(res$xopt)$deltahat
InvDeltahat <- solve(Deltahat)
temp0 <- -(n*p/2)*log(2*pi);
temp1 <- -(n/2)*log(det(Deltahat)); 
temp2 <- -(n/2)*Trace(InvDeltahat%*%Sres);
numdir = 2
if (numdir < p) temp3 <- -(n/2)*sum(eigen(InvDeltahat%*%Sres)$values[(numdir+1):p]);
loglik <- Re(temp0 + temp1 + temp2 + temp3);
print(-loglik)

# [1] 4062.263
