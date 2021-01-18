#' Manifold optimization
#'
#' Optimize a function on a manifold.
#'
#' @param prob \link{Problem definition}
#' @param mani.defn Either a \link{Product manifold definition} or one of the
#'        \link{Manifold definitions}
#' @param method Name of optimization method. Currently supported methods are:
#' \itemize{
#' \item{\code{"LRBFGS"}}: Limited-memory RBFGS
#' \item{\code{"LRTRSR1"}}: Limited-memory RTRSR1
#' \item{\code{"RBFGS"}}: Riemannian BFGS
#' \item{\code{"RBroydenFamily"}}: Riemannian Broyden family
#' \item{\code{"RCG"}}: Riemannian conjugate gradients
#' \item{\code{"RNewton"}}: Riemannian line-search Newton
#' \item{\code{"RSD"}}: Riemannian steepest descent
#' \item{\code{"RTRNewton"}}: Riemannian trust-region Newton
#' \item{\code{"RTRSD"}}: Riemannian trust-region steepest descent
#' \item{\code{"RTRSR1"}}: Riemannian trust-region symmetric rank-one update
#' \item{\code{"RWRBFGS"}}: Riemannian BFGS
#' }
#' See Huang et al (2016a, 2016b) for details.
#' @param mani.params Arguments to configure the manifold. Construct with
#'        \link{get.manifold.params}
#' @param solver.params Arguments to configure the solver. Construct with
#'        \link{get.solver.params}
#' @param deriv.params Arguments to configure numerical differentiation for
#'        gradient and Hessian, which are used if those functions are not
#'        specified. Construct with \link{get.deriv.params}
#' @param x0 Starting point for optimization. A numeric vector whose dimension
#'        matches the total dimension of the overall problem
#' @param H0 Initial value of Hessian. A \eqn{d \times d} matrix, where \eqn{d}
#'        is the dimension of \code{x0}
#' @param has.hhr Indicates whether to apply the idea in Huang et al
#'        (2015) section 4.1 (TRUE or FALSE)
#'
#' @return
#' \item{xopt}{Point returned by the solver}
#' \item{fval}{Value of the function evaluated at \code{xopt}}
#' \item{normgf}{Norm of the final gradient}
#' \item{normgfgf0}{Norm of the gradient at final iterate divided by norm of
#'   the gradient at initiate iterate}
#' \item{iter}{Number of iterations taken by the solver}
#' \item{num.obj.eval}{Number of function evaluations}
#' \item{num.grad.eval}{Number of gradient evaluations}
#' \item{nR}{Number of retraction evaluations}
#' \item{nV}{Number of occasions in which vector transport is first computed}
#' \item{nVp}{Number of remaining computations of vector transport (excluding count in nV)}
#' \item{nH}{Number of actions of Hessian}
#' \item{elapsed}{Elapsed time for the solver (in seconds)}
#'
#' @details \code{moptim} is an alias for  \code{manifold.optim}.
#'
#' @examples
#' @example inst/examples/brockett/rproblem/driver-minimal.Rin
#' @example inst/examples/brockett/cpp_sourceCpp/driver-minimal.Rin
#'
#' @references
#'
#' Wen Huang, P.A. Absil, K.A. Gallivan, Paul Hand (2016a). "ROPTLIB: an
#' object-oriented C++ library for optimization on Riemannian manifolds."
#' Technical Report FSU16-14, Florida State University.
#'
#' Wen Huang, Kyle A. Gallivan, and P.A. Absil (2016b).
#' Riemannian Manifold Optimization Library.
#' URL \url{http://www.math.fsu.edu/~whuang2/pdf/USER_MANUAL_for_2016-04-29.pdf}
#'
#' Wen Huang, K.A. Gallivan, and P.A. Absil (2015). A Broyden Class of
#' Quasi-Newton Methods for Riemannian Optimization. SIAM  Journal on
#' Optimization, 25(3):1660-1685.
#'
#' S. Martin, A. Raim, W. Huang, and K. Adragni (2020). "ManifoldOptim: 
#' An R Interface to the ROPTLIB Library for Riemannian Manifold Optimization."
#' Journal of Statistical Software, 93(1):1-32.
#'
#' @name manifold.optim
manifold.optim <- function(prob, mani.defn, method = "LRBFGS", 
	mani.params = get.manifold.params(), solver.params = get.solver.params(),
	deriv.params = get.deriv.params(), x0 = NULL, H0 = NULL, has.hhr = FALSE)
{
	if (class(mani.params) != "ManifoldOptim::manifold.params") {
		stop("mani.params argument should be of class ManifoldOptim::manifold.params")
	}
	
	if (class(solver.params) != "ManifoldOptim::solver.params") {
		stop("solver.params argument should be of class ManifoldOptim::solver.params")
	}

	if (class(deriv.params) != "ManifoldOptim::deriv.params") {
		stop("deriv.params argument should be of class ManifoldOptim::deriv.params")
	}

	if (is.null(x0)) {
		x0 <- matrix(0, 0, 0)
	}

	if (is.null(H0)) {
		H0 <- matrix(0, 0, 0)
	}

	if (is.null(solver.params$IsCheckParams)) { solver.params$IsCheckParams = 0 }
	if (is.null(solver.params$IsCheckGradHess)) { solver.params$IsCheckGradHess = 0 }
	if (is.null(mani.params$IsCheckParams)) { mani.params$IsCheckParams = 0 }

	if (class(mani.defn) == "ManifoldOptim::manifold.defn") {
		mani.defn.list <- get.product.defn(mani.defn)
	} else if (class(mani.defn) == "ManifoldOptim::product.manifold.defn") {
		mani.defn.list <- mani.defn
	} else {
		stop("mani.defn is not a valid manifold definition")
	}

	res <- .Call('manifold_optim', PACKAGE = 'ManifoldOptim', x0, H0, prob,
		mani.defn.list, mani.params, solver.params, deriv.params, method, has.hhr)
	class(res) <- "ManifoldOptim"
	return(res)
}

#' @name manifold.optim
moptim <- manifold.optim

