#' Get parameters to initialize solver
#'
#' @details
#' Solver-specific parameters may also be added to the object returned
#' from \code{get.solver.params}, via standard list manipulation. Interested
#' users should refer to Huang et al (2016b) for available options.
#'
#' @param isconvex Indicator for whether the function is convex (TRUE or FALSE)
#' @param DEBUG Verbosity level in \{0,1,2,3\}. Use 0 for quietest with no
#'        messages printed. Use 3 for most verbose,
#' @param Tolerance Tolerance used to assess convergence. See Huang et al
#'        (2016b) for details on how this is used,
#' @param Max_Iteration Maximum iterations to be used by the solver
#'        (a non-negative integer),
#' @param IsCheckParams Should solver check inputs and print summary message
#'        before optimization (TRUE or FALSE),
#' @param IsCheckGradHess Check correctness of the gradient and Hessian
#'        functions (TRUE or FALSE).
#' @param ... Additional arguments to pass to the solver. These are not
#'        validated by the \code{get.solver.params} function. Users should 
#'        refer to the C++ library's user manual for available arguments.
#'
#' @return List containing input arguments for solver
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
#'
get.solver.params <- function(isconvex = FALSE, DEBUG = 0, Tolerance = 1e-4,
	Max_Iteration = 1000, IsCheckParams = FALSE, IsCheckGradHess = FALSE, ...) {
	if (!is.logical(isconvex)) { stop("isconvex must be logical") }
	if (DEBUG %notin% c(0,1,2,3)) { stop("DEBUG must be one of {0,1,2,3}") }
	if (!is.numeric(Tolerance)) { stop("Tolerance must be numeric") }
	if (!is.numeric(Max_Iteration) || Max_Iteration <= 0) { stop("Max_Iteration must be a positive integer") }
	if (!is.logical(IsCheckParams)) { stop("IsCheckParams must be logical") }
	if (!is.logical(IsCheckGradHess)) { stop("IsCheckGradHess must be logical") }
	structure(
		list(isconvex = as.integer(isconvex),
			 DEBUG = as.integer(DEBUG),
			 Tolerance = as.numeric(Tolerance),
			 Max_Iteration = as.integer(Max_Iteration),
			 IsCheckParams = as.integer(IsCheckParams),
			 IsCheckGradHess = as.integer(IsCheckGradHess),
			 ...),
		class = c("ManifoldOptim::solver.params")
	)
}

#' Get parameters to initialize manifold
#'
#' @param IsCheckParams Should internal manifold object check inputs and
#'        print summary message before optimization (TRUE or FALSE)
#'
#' @return List containing input arguments for manifold
#'
get.manifold.params <- function(IsCheckParams = FALSE) {
	if (!is.logical(IsCheckParams)) { stop("IsCheckParams must be logical") }
	structure(
		list(IsCheckParams = as.integer(IsCheckParams)),
		class = c("ManifoldOptim::manifold.params")
	)
}

#' Get parameters to initialize numerical differentiation
#'
#' @param EpsNumericalGrad The "epsilon" used to perturb the objective functon
#'        when computing numerical gradients
#' @param EpsNumericalHessEta The "epsilon" used to perturb the objective functon
#'        when computing numerical HessEta
#'
#' @return List containing input arguments for numerical differentiation
#'
get.deriv.params <- function(EpsNumericalGrad = 1e-6, EpsNumericalHessEta = 1e-4) {
	if (!is.numeric(EpsNumericalGrad)) { stop("EpsNumericalGrad must be numeric") }
	if (!is.numeric(EpsNumericalHessEta)) { stop("EpsNumericalHessEta must be numeric") }
	structure(
		list(EpsNumericalGrad = as.numeric(EpsNumericalGrad),
			EpsNumericalHessEta = as.numeric(EpsNumericalHessEta)),
		class = c("ManifoldOptim::deriv.params")
	)
}

