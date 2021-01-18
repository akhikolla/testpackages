# Date Created: 5/10/16
# Description:
# This is a simulation of an envelope method applied to logistic regression as 
# detailed on page 23 of:
# D. Cook, X. Zhang "Foundations for Envelope Models and Methods." 2014.

require(MASS)
require(ldr)
library(Rcpp)
library(RcppArmadillo)
library(ManifoldOptim)
library(mvtnorm)

#set.seed(1234)  # seed 1
set.seed(7234)   # seed 2

# number of data points to be simulated
#n <- 300 
n <- 150

# true beta
beta <- matrix(c(0.25, 0.25), nr=2);

# normalized beta for first eigenvector
nbeta <- beta / norm(beta, "F")

# find complementary eigenvector to form covariance
nbeta_comp <- diag(2) - nbeta %*% t(nbeta)
evals <- eigen(nbeta_comp)
e_vec <- matrix(evals$vectors[,1], nrow = 2)

# construct a basis
sig_basis <- cbind(nbeta,e_vec)

# paper defined eigenvalues
paper_evals <- matrix(diag(c(10, .1)), nr=2);

# construct a covariance
Sigma <- sig_basis %*% paper_evals %*% t(sig_basis)

# construct predictors
X <- t(mvrnorm(n=n, mu=c(0,0), Sigma=Sigma));

# construct responses
Y <- c()

for(i in 1:n){
  
  # construct the logit value
  alpha <- t(beta) %*% X[,i];
  la <- 1 / (1 + exp(-alpha)) 
  #print(la)

  # generate data from Bernoulli  
  #Y <- rbinom(150, 1, pr);
  Y <- c(Y, rbinom(1, 1, la));

}

# construct a matrix from y
Y <- matrix(Y, nr=n);

# fit a logistic regression
m1 <- glm(as.factor(Y) ~ t(X), family="binomial")
summary(m1)

# look at estimates of beta
m1$coefficients

# setup dimensions (note: r = p, this is a difference between papers)
r <- 2
u <- 1

# marginal covariance of predictors
Sx <- cov(t(X));

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
	
	     # back out Xi
	     xi <- matrix(X[,i], nr=2);

	     # and back out Yi
	     yi <- Y[i];

	     # estimate beta
	     betahat <- gamma_hat %*% eta_hat;

	     # Calculate vi from page 14
	     vi <- alpha_hat + t(betahat) %*% xi;

	     # Calculate A(vi) from Table 1 on page 14
	     A_vi <- 1 + exp(vi);

	     # Calculate C(vi) from Table 1 on page 14
	     C_vi <- yi %*% vi - log(A_vi);

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

# At this point do hill climbing with random restarts to avoid 
# convergence to local minima
min_fval <- Inf;
min_xopt <- c();

num_restarts <- 20
for(restart in 1:num_restarts){

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
      solver.params <- get.solver.params(DEBUG = 2, Tolerance = 1e-300, Max_Iteration = 1000, IsCheckParams = TRUE, IsCheckGradHess = FALSE)

      mani.defn <- get.product.defn(get.grassmann.defn(r,u), get.euclidean.defn(u,1), get.euclidean.defn(u,1))

      res <- manifold.optim(prob, mani.defn, method = "LRBFGS", mani.params = mani.params, solver.params = solver.params, x0 = x0)


      # look at the output and save off if this is the best so far
      print(res)

      if(res$fval < min_fval){
         min_fval <- res$fval
	 min_xopt <- res$xopt
      }

}


#method = "RTRSR1", 

par <- tx(min_xopt) 
print(par)
print(min_fval)

cat("\nactual beta:\n")
print(beta)

cat("\nestimate of beta from glm:\n")
print(m1$coefficients[2:3])

cat("\nestimate of beta from envelope:\n")
print(par$gammahat %*% par$etahat)

# SRM NOTE: 
# using LRBFGS algorithm on seed 2, final obj value is 111.437 and betahat_env
# is:
#> print(par$gammahat %*% par$etahat)
#          [,1]
#[1,] 0.2463228
#[2,] 0.2405503

# same setup on seed 1 with n=300:
#> print(par$gammahat %*% par$etahat)
#          [,1]
#[1,] 0.2521251
#[2,] 0.2531871


