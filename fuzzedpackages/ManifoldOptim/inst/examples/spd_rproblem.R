library(ManifoldOptim)
library(mvtnorm)

set.seed(1234)

n <- 2000
p <- 5
mu.true <- rep(0,p)
Sigma.true <- diag(2,p) + 0.01
y <- rmvnorm(n, mean = mu.true, sigma = Sigma.true)

# It looks like optimizer can return matrices that are not quite symmetric.
# dmvnorm does not accept this, so let's symmetrize before we evaluate.

tx <- function(x) {
	S <- matrix(x, p, p)
	S[lower.tri(S)] <- t(S)[lower.tri(S)]
	return(S)
}
f <- function(x) {
	Sigma <- tx(x)
	-sum(dmvnorm(y, mean = mu.true, sigma = Sigma, log = TRUE))
}

# Take the problem above (defined in R) and make a Problem object for it.
# The Problem can be passed to the optimizer in C++
mod <- Module("ManifoldOptim_module", PACKAGE = "ManifoldOptim")
prob <- new(mod$RProblem, f)

X0 <- diag(1, p)
x0 <- as.numeric(X0)

# Test the obj and grad fn
prob$objFun(x0)
# head(prob$gradFun(x0))

mani.params <- get.manifold.params(IsCheckParams = TRUE)
solver.params <- get.solver.params(DEBUG = 0, Tolerance = 1e-4,
	Max_Iteration = 1000, IsCheckParams = TRUE, IsCheckGradHess = FALSE)
mani.defn <- get.spd.defn(p)

res <- manifold.optim(prob, mani.defn, method = "LRBFGS", 
	mani.params = mani.params, solver.params = solver.params, x0 = x0)
print(res)
head(tx(res$xopt))

