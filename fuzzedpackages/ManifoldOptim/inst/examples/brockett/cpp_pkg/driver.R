library(ManifoldOptim)

set.seed(1234)

p <- 5
n <- 150
B <- matrix(rnorm(n*n), nrow=n)
B <- B + t(B) # force symmetric
D <- diag(p:1, p)

tx <- function(x) { matrix(x, n, p) }

# The Problem class is written in C++. Get a handle to it and set it up from R
mod <- Module("Brockett_module", PACKAGE = "ManifoldOptim")
prob <- new(mod$BrockettProblem, B, D)

X0 <- orthonorm(matrix(rnorm(n*p), nrow=n, ncol=p))
x0 <- as.numeric(X0)
eta <- diag(1, n*p, 1)
prob$objFun(x0)						# Test the obj fn
head(tx(prob$gradFun(x0)))			# Test the grad fn
# head(prob$hessEtaFun(x0, eta))	# Test the Hess fn (numerical deriv is very slow)
prob$GetB()[1:5,1:5]
head(prob$GetD())

# ----- Run manifold.optim -----
mani.params <- get.manifold.params(IsCheckParams = TRUE)
solver.params <- get.solver.params(DEBUG = 0, Tolerance = 1e-4,
	Max_Iteration = 1000, IsCheckParams = TRUE, IsCheckGradHess = FALSE)
mani.defn <- get.stiefel.defn(n, p)

res <- manifold.optim(prob, mani.defn, method = "RTRSR1",
	mani.params = mani.params, solver.params = solver.params, x0 = x0)
print(res)
head(tx(res$xopt))
