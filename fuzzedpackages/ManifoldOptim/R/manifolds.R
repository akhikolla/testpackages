#' Product manifold definition
#'
#' Define a product manifold composed of simpler manifolds
#'
#' @param ... One or more simpler \link{Manifold definitions}
#'
#' @return List containing manifold definitions for the product manifold
#'
#' @name Product manifold definition
#'
#' @examples 
#' mani.defn1 <- get.product.defn(get.sphere.defn(n=5), get.spd.defn(n=5))
#' mani.defn2 <- get.product.defn(
#'     get.stiefel.defn(n=10, p=5),
#'     get.stiefel.defn(n=7, p=3),
#'     get.grassmann.defn(n=10, p=5)
#' )
#' 
#' \dontrun{
#' # --- Estimate jointly: Sigma in SPD manifold and mu in sphere manifold ---
#' library(mvtnorm)
#' n <- 400
#' p <- 3
#' mu.true <- rep(1/sqrt(p), p)
#' Sigma.true <- diag(2,p) + 0.1
#' y <- rmvnorm(n, mean = mu.true, sigma = Sigma.true)
#' 
#' tx <- function(x) {
#'     idx.mu <- 1:p
#'     idx.S <- 1:p^2 + p
#'     mu <- x[idx.mu]
#'     S <- matrix(x[idx.S], p, p)
#'     list(mu = mu, Sigma = S)
#' }
#' f <- function(x) {
#'     par <- tx(x)
#'     -sum(dmvnorm(y, mean = par$mu, sigma = par$Sigma, log = TRUE))
#' }
#' 
#' mod <- Module("ManifoldOptim_module", PACKAGE = "ManifoldOptim")
#' prob <- new(mod$RProblem, f)
#' 
#' mu0 <- diag(1, p)[,1]
#' Sigma0 <- diag(1, p)
#' x0 <- c(mu0, as.numeric(Sigma0))
#' 
#' mani.defn <- get.product.defn(get.sphere.defn(p), get.spd.defn(p))
#' mani.params <- get.manifold.params()
#' solver.params <- get.solver.params(isconvex = TRUE)
#' 
#' res <- manifold.optim(prob, mani.defn, method = "LRBFGS", 
#'     mani.params = mani.params, solver.params = solver.params, x0 = x0)
#' }
#'
get.product.defn <- function(...) {
	x <- list(...)
	L <- length(x)
	if (L == 0) { stop("Must provide at least one simple manifold") }
	
	for (i in 1:L) {
		if (class(x[[i]]) != "ManifoldOptim::manifold.defn") {
			stop("Each entry of mani.defn argument should be of class ManifoldOptim::manifold.defn")
		}
	}
	
	structure(x, class = c("ManifoldOptim::product.manifold.defn"))
}

#' Manifold definitions
#'
#' Get definitions for simple manifolds
#'
#' @details The functions define manifolds as follows:
#' \itemize{
#' \item{\code{get.stiefel.defn}}: Stiefel manifold
#'    \eqn{\{X \in R^{n \times p} : X^T X = I\}}
#' \item{\code{get.grassmann.defn}}: Grassmann manifold of \eqn{p}-dimensional
#'   subspaces in \eqn{R^n}
#' \item{\code{get.spd.defn}}: Manifold of \eqn{n \times n} symmetric positive
#'    definite matrices
#' \item{\code{get.sphere.defn}}: Manifold of \eqn{n}-dimensional vectors on
#'    the unit sphere
#' \item{\code{get.euclidean.defn}}: Euclidean \eqn{R^{n \times m}} space
#' \item{\code{get.lowrank.defn}}: Low-rank manifold
#'    \eqn{\{ X \in R^{n \times m} : \textrm{rank}(X) = p \}}
#' \item{\code{get.orthgroup.defn}}: Orthonormal group
#'    \eqn{\{X \in R^{n \times n} : X^T X = I\}}
#' }
#'
#' @param n Dimension for manifold object (see Details)
#' @param p Dimension for manifold object (see Details)
#' @param m Dimension for manifold object (see Details)
#' @param numofmani Multiplicity of this space. For example, use
#'        \code{numofmani = 2} if problem requires 2 points from this manifold
#' @param ParamSet A positive integer indicating a set of properties for the
#'        manifold which can be used by the solver. See Huang et al (2016b)
#'        for details.
#'
#' @return List containing input arguments and name field denoting the type of manifold
#'
#' @name Manifold definitions
#'
#' @references
#' Wen Huang, P.A. Absil, K.A. Gallivan, Paul Hand (2016a). "ROPTLIB: an
#' object-oriented C++ library for optimization on Riemannian manifolds."
#' Technical Report FSU16-14, Florida State University.
#'
#' Wen Huang, Kyle A. Gallivan, and P.A. Absil (2016b).
#' Riemannian Manifold Optimization Library.
#' URL \url{http://www.math.fsu.edu/~whuang2/pdf/USER_MANUAL_for_2016-04-29.pdf}
#'
#' S. Martin, A. Raim, W. Huang, and K. Adragni (2020). "ManifoldOptim: 
#' An R Interface to the ROPTLIB Library for Riemannian Manifold Optimization."
#' Journal of Statistical Software, 93(1):1-32.
NULL


#' @name Manifold definitions
get.stiefel.defn <- function(n, p, numofmani = 1L, ParamSet = 1L) {
	if (!is.numeric(n)) { stop("n must be a number") }
	if (!is.numeric(p)) { stop("p must be a number") }
	if (!is.numeric(numofmani)) { stop("numofmani must be a number") }
	if (!is.numeric(ParamSet)) { stop("ParamSet must be a number") }
	structure(
		list(name = "Stiefel",
			 n = as.integer(n),
			 p = as.integer(p),
			 m = 1L,
			 numofmani = as.integer(numofmani),
			 ParamSet = as.integer(ParamSet)),
		class = c("ManifoldOptim::manifold.defn")
	)
}

#' @name Manifold definitions
get.grassmann.defn <- function(n, p, numofmani = 1L, ParamSet = 1L) {
	if (!is.numeric(n)) { stop("n must be a number") }
	if (!is.numeric(p)) { stop("p must be a number") }
	if (!is.numeric(numofmani)) { stop("numofmani must be a number") }
	if (!is.numeric(ParamSet)) { stop("ParamSet must be a number") }
	structure(
		list(name = "Grassmann",
			 n = as.integer(n),
			 p = as.integer(p),
			 m = 1L,
			 numofmani = as.integer(numofmani),
			 ParamSet = as.integer(ParamSet)),
		class = c("ManifoldOptim::manifold.defn")
	)
}

#' @name Manifold definitions
get.spd.defn <- function(n, numofmani = 1L, ParamSet = 1L) {
	if (!is.numeric(n)) { stop("n must be a number") }
	if (!is.numeric(numofmani)) { stop("numofmani must be a number") }
	if (!is.numeric(ParamSet)) { stop("ParamSet must be a number") }
	structure(
		list(name = "SPDManifold",
			 n = as.integer(n),
			 p = as.integer(n),
			 m = 1L,
			 numofmani = as.integer(numofmani),
			 ParamSet = as.integer(ParamSet)),
		class = c("ManifoldOptim::manifold.defn")
	)
}

#' @name Manifold definitions
get.sphere.defn <- function(n, numofmani = 1L, ParamSet = 1L) {
	if (!is.numeric(n)) { stop("n must be a number") }
	if (!is.numeric(numofmani)) { stop("numofmani must be a number") }
	if (!is.numeric(ParamSet)) { stop("ParamSet must be a number") }
	structure(
		list(name = "Sphere",
			 n = as.integer(n),
			 p = as.integer(n),
			 m = 1L,
			 numofmani = as.integer(numofmani),
			 ParamSet = as.integer(ParamSet)),
		class = c("ManifoldOptim::manifold.defn")
	)
}

# @name Manifold definitions
#get.oblique.defn <- function(n, m, numofmani = 1L, ParamSet = 1L) {
# if (!is.numeric(n)) { stop("n must be a number") }
# if (!is.numeric(m)) { stop("m must be a number") }
# if (!is.numeric(numofmani)) { stop("numofmani must be a number") }
# if (!is.numeric(ParamSet)) { stop("ParamSet must be a number") }
#	structure(
#		list(name = "Oblique",
#			 n = as.integer(n),
#			 p = 1L,
#			 m = as.integer(m),
#			 numofmani = as.integer(numofmani),
#			 ParamSet = as.integer(ParamSet)),
#		class = c("ManifoldOptim::manifold.defn")
#	)
#}

#' @name Manifold definitions
get.euclidean.defn <- function(n, m, numofmani = 1L, ParamSet = 1L) {
	if (!is.numeric(n)) { stop("n must be a number") }
	if (!is.numeric(m)) { stop("m must be a number") }
	if (!is.numeric(numofmani)) { stop("numofmani must be a number") }
	if (!is.numeric(ParamSet)) { stop("ParamSet must be a number") }
	structure(
		list(name = "Euclidean",
			 n = as.integer(n),
			 p = 1L,
			 m = as.integer(m),
			 numofmani = as.integer(numofmani),
			 ParamSet = as.integer(ParamSet)),
		class = c("ManifoldOptim::manifold.defn")
	)
}

#' @name Manifold definitions
get.lowrank.defn <- function(n, m, p, numofmani = 1L, ParamSet = 1L) {
	if (!is.numeric(n)) { stop("n must be a number") }
	if (!is.numeric(m)) { stop("m must be a number") }
	if (!is.numeric(p)) { stop("p must be a number") }
	if (!is.numeric(numofmani)) { stop("numofmani must be a number") }
	if (!is.numeric(ParamSet)) { stop("ParamSet must be a number") }
	structure(
		list(name = "LowRank",
			 n = as.integer(n),
			 p = as.integer(p),
			 m = as.integer(m),
			 numofmani = as.integer(numofmani),
			 ParamSet = as.integer(ParamSet)),
		class = c("ManifoldOptim::manifold.defn")
	)
}

#' @name Manifold definitions
get.orthgroup.defn <- function(n, numofmani = 1L, ParamSet = 1L) {
	if (!is.numeric(n)) { stop("n must be a number") }
	if (!is.numeric(numofmani)) { stop("numofmani must be a number") }
	if (!is.numeric(ParamSet)) { stop("ParamSet must be a number") }
	structure(
		list(name = "OrthGroup",
			 n = as.integer(n),
			 p = 1L,
			 m = 1L,
			 numofmani = as.integer(numofmani),
			 ParamSet = as.integer(ParamSet)),
		class = c("ManifoldOptim::manifold.defn")
	)
}

# @name Manifold definitions
#get.l2sphere.defn <- function(n, numofmani = 1L, ParamSet = 1L) {
#if (!is.numeric(n)) { stop("n must be a number") }
#if (!is.numeric(numofmani)) { stop("numofmani must be a number") }
#if (!is.numeric(ParamSet)) { stop("ParamSet must be a number") }
#	structure(
#		list(name = "L2Sphere",
#			 n = as.integer(n),
#			 p = 1L,
#			 m = 1L,
#			 numofmani = as.integer(numofmani),
#			 ParamSet = as.integer(ParamSet)),
#		class = c("ManifoldOptim::manifold.defn")
#	)
#}

