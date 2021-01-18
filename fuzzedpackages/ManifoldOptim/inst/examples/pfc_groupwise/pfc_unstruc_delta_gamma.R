require(MASS)
require(ldr)
library(Rcpp)
library(RcppArmadillo)
library(ManifoldOptim)
library(mvtnorm)
source("MainFunctions_P.R")

set.seed(1234)

# ---- Start data generation ----

# start by constructing an unstructured covariance amoung predictors 
n = 200; p = 12;
y <- array( rnorm(n=n, sd=3), c(n, 1));
fy <- scale(cbind(y, abs(y)), TRUE, FALSE)
A = matrix(rnorm(p^2), nc=p); D = t(A)%*%A
Gamma = eigen(D)$vectors[,1:2]

# variances for all predictors
sigma2 <- 1

# build true Delta
Delta = matrix(rnorm(p^2),nc=p); 
Delta = t(Delta)%*%Delta;
ev = eigen(Delta);
Delta = t(ev$vectors) %*% diag(seq(1,p)) %*% ev$vectors

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

	temp2 <- -(n/2)*Trace(InvDeltahat %*% Sres) - (n/2)*Trace(Gammahat)

	# check for singular matrices
	#tmp_test <- t(Gammahat) %*% InvDeltahat %*% Gammahat
	#print(tmp_test)


	# Kofi's equation
	#temp2 <- -(n/2)*Trace(Sigmahat %*% InvDeltahat) +(n/2)*Trace(t(Gammahat) %*% InvDeltahat %*% Sfit %*% Gammahat %*% solve(diag(1e-12,2) + t(Gammahat) %*% InvDeltahat %*% Gammahat))
	

	loglik <- Re(temp0 + temp1 + temp2);
	return(-loglik)

}

# additional metric for looking at angular distances
dist.space <- function(A, B)
{
  p <- NROW(A)
  stopifnot(p == NROW(B))
  d <- NCOL(A)
  D <- NCOL(B)
  II <- diag(1,p)
  
  if (d < D) return(Matrix::norm( (II - B %*% t(B)) %*% A, type = "F"))  
  else return(Matrix::norm( (II - A %*% t(A)) %*% B, type = "F"))
}




# Take the problem above (defined in R) and make a Problem object for it.
# The Problem can be passed to the optimizer in C++
mod <- Module("ManifoldOptim_module", PACKAGE = "ManifoldOptim")
prob <- new(mod$RProblem, f)

# initialize with PFC output
#x0 = c(as.numeric(fit0$Deltahat), as.numeric(fit0$Gammahat))

# Test the obj and grad fn
#prob$objFun(x0)
#head(prob$gradFun(x0))

# At this point do hill climbing with random restarts to avoid 
# convergence to local minima
min_fval <- Inf;
min_xopt <- c();

num_restarts <- 10
for(restart in 1:num_restarts){



# initialize randomly 
A = matrix(rnorm(p^2), nc=p); 
D = t(A)%*%A;
evals = eigen(D);
gamma0 <- matrix(evals$vectors[,1:2],nr=p);
delta0 = matrix(rnorm(p^2),nc=p); 
delta0 = t(delta0)%*%delta0;
ev = eigen(delta0);
delta0 <- t(ev$vectors) %*% diag(seq(1,p)) %*% ev$vectors
x0 = c(as.numeric(delta0), as.numeric(gamma0))


mani.params <- get.manifold.params(IsCheckParams = TRUE)
solver.params <- get.solver.params(DEBUG = 2, Tolerance = 1e-12,
	Max_Iteration = 1000, IsCheckParams = TRUE, IsCheckGradHess = FALSE)
mani.defn <- get.product.defn(get.spd.defn(p), get.grassmann.defn(p,2))

res <- manifold.optim(prob, mani.defn, method = "LRBFGS", 
	mani.params = mani.params, solver.params = solver.params, 
	x0 = x0)

print(res)

      if(res$fval < min_fval){
         min_fval <- res$fval
	 min_xopt <- res$xopt
      }

}



#print(tx(res$xopt))


# examine likelihood
prob$objFun(res$xopt)

#[1] 2949.486



# examine heatmap
heatmap(tx(res$xopt)$deltahat, Rowv=NA, Colv=NA)

# examine sum of eigen values in diff
Trace(fit0$Deltahat - Delta)
#[1] 0.3662055
Trace(tx(min_xopt)$deltahat - Delta)
#[1] 0.4785665

# examine the sum of singular values in diff
sum(svd(Gamma - fit0$Gammahat)$d)
#[1] 3.994087
sum(svd(Gamma - tx(min_xopt)$gammahat)$d)
#[1] 3.07502

dist.space(fit0$Deltahat, Delta)
dist.space(tx(min_xopt)$deltahat, Delta)

dist.space(fit0$Gammahat, Gamma)
dist.space(tx(min_xopt)$gammahat, Gamma)

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
