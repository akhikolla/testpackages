library(ManifoldOptim)

set.seed(1234)

p <- 25
d <- 2
a <- matrix(rnorm(p^2), ncol = p)
A <- t(a) %*% a

tx <- function(x) {
	matrix(x, p, d)
}
f <- function(x) {
	X <- tx(x)
	Y <- X[,1:d]
	0.5 * sum(diag(t(Y) %*% A %*% Y))
}
g <- function(x) {
	X <- tx(x)
	t(X) %*% A %*% X
}

# Create an RProblem for the above
mod <- Module("ManifoldOptim_module", PACKAGE = "ManifoldOptim")
prob <- new(mod$RProblem, f, g)

# Random starting matrix
m <- matrix(rnorm(p^2), ncol=p)
M <- t(m) %*% m
W <- list(Qt=eigen(M)$vectors, dim=c(d,p), A=A)
decomp <- svd(W$Qt, nu=2, nv=2)
X0 <- decomp$u
x0 <- as.numeric(X0)

# Test the obj and grad fn
prob$objFun(x0)
head(prob$gradFun(x0))

mani.params <- get.manifold.params(IsCheckParams = TRUE)
solver.params <- get.solver.params(DEBUG = 0, Tolerance = 1e-4,
	Max_Iteration = 1000, IsCheckParams = TRUE, IsCheckGradHess = TRUE)
mani.defn <- get.grassmann.defn(p, d)

# TBD: R crashes if dimensions don't match up

res <- manifold.optim(prob, mani.defn, method = "LRBFGS", 
	mani.params = mani.params, solver.params = solver.params, x0 = x0)
print(res)
head(tx(res$xopt))

