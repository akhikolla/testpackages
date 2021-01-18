library(ManifoldOptim)
library(mvtnorm)

set.seed(1234)

n <- 2000
p <- 10
mu.true <- rep(1/sqrt(p), p)
Sigma.true <- diag(2,p) + 0.01
y <- rmvnorm(n, mean = mu.true, sigma = Sigma.true)

f <- function(x) {
	mu <- x
	-sum(dmvnorm(y, mean = mu, sigma = Sigma.true, log = TRUE))
}

# Take the problem above (defined in R) and make a Problem object for it.
# The Problem can be passed to the optimizer in C++
mod <- Module("ManifoldOptim_module", PACKAGE = "ManifoldOptim")
prob <- new(mod$RProblem, f)

x0 <- diag(1, p)[,1]

# Test the obj and grad fn
prob$objFun(x0)
# head(prob$gradFun(x0))

mani.params <- get.manifold.params(IsCheckParams = TRUE)
solver.params <- get.solver.params(DEBUG = 0, Tolerance = 1e-4,
	Max_Iteration = 1000, IsCheckParams = TRUE, IsCheckGradHess = FALSE)
mani.defn <- get.sphere.defn(p)

res <- manifold.optim(prob, mani.defn, method = "LRBFGS", 
	mani.params = mani.params, solver.params = solver.params, x0 = x0)
print(res)
print(res$xopt)
