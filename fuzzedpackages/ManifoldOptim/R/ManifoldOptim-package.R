#' @title
#' Problem definition
#'
#' @description
#' Define a problem for ManifoldOptim to solve.
#'
#' @details
#' A problem definition contains an objective function \eqn{f} and a gradient
#' function \eqn{g}. The gradient \eqn{g} is computed as if \eqn{f} is defined
#' on a Euclidean space. If \eqn{g} is not specified it will be computed
#' numerically, which is potentially much slower.
#'
#' The easiest way to define a problem is completely in R. Example 1
#' below illustrates how to construct a problem using a given \eqn{f} and
#' \eqn{g}. Example 2 constructs the same problem without providing \eqn{g}.
#' The \code{Rcpp Module} framework (Eddelbuettel, 2013) creates underlying
#' C++ objects necessary to invoke the \code{ROPTLIB} library.
#'
#' The performance of solving an \code{RProblem} may be too slow for some
#' applications; here, the C++ optimizer calls R functions,
#' which requires some overhead. A faster alternative is to code your problem
#' in C++ directly, and allow it to be manipulated in R. An
#' example is provided in this package, under
#' \code{tests/brockett/cpp_standalone/}. Example 3 below shows how to
#' instantiate this problem.
#'
#' Package authors may want to use \code{ManifoldOptim} within a package to solve
#' a problem written in C++. In this case, the author would probably
#' not want to use \code{sourceCpp}, but instead have the problem compiled
#' when the package was installed. An example is provided within this package;
#' \code{tests/brockett/cpp_pkg/driver.R} instantiates the problem defined in:
#'
#' \code{src/ManifoldOptim/BrockettProblem.cpp}.
#'
#' @examples
#' \dontrun{
#' # --- Example 1: Define a problem in R ---
#' f <- function(x) { ... }
#' g <- function(x) { ... }
#' mod <- Module("ManifoldOptim_module", PACKAGE = "ManifoldOptim")
#' prob <- new(mod$RProblem, f, g)
#'
#' # --- Example 2: Define a problem in R without specifying gradient ---
#' f <- function(x) { ... }
#' mod <- Module("ManifoldOptim_module", PACKAGE = "ManifoldOptim")
#' prob <- new(mod$RProblem, f)
#'
#' # --- Example 3: Instantiate a problem written in C++ ---
#' p <- 5; n <- 150
#' B <- matrix(rnorm(n*n), nrow=n)
#' B <- B + t(B) # force symmetric
#' D <- diag(p:1, p)
#' Rcpp::sourceCpp("brockett_problem.cpp")
#' prob <- new(BrockettProblem, B, D)
#' }
#'
#' @name Problem definition
#'
#' @references
#' Dirk Eddelbuettel. Seamless R and C++ Integration with Rcpp,
#'   Chapter 7: Modules, pages 83-102. Springer New York, New York, NY, 2013.
#'
#' Wen Huang, P.A. Absil, K.A. Gallivan, Paul Hand (2016a). "ROPTLIB: an
#' object-oriented C++ library for optimization on Riemannian manifolds."
#' Technical Report FSU16-14, Florida State University.
#'
#' S. Martin, A. Raim, W. Huang, and K. Adragni (2020). "ManifoldOptim: 
#' An R Interface to the ROPTLIB Library for Riemannian Manifold Optimization."
#' Journal of Statistical Software, 93(1):1-32.
NULL

#' @title
#' Overview of important files.
#' 
#' @description
#' Internal design of the ManifoldOptim portion of the embedded C++ code.
#' Most ManifoldOptim users should not need this.
#' ROPTLIB source code is also included in this package, but is not described here;
#' see Huang et al (2016a) for documentation on that portion of the code.
#' 
#' @details
#' \describe{
#' \item{\code{src/ManifoldOptim/BrockettProblem.cpp}:}{
#' The Brockett problem, written as a module that can be invoked from within
#' the ManifoldOptim package. This serves as an example for package authors
#' who wish to expose modules to their users. Code to invoke this example from
#' outside of the ManifoldOptim package is provided in 
#' \code{inst/examples/brockett/cpp_pkg}.
#' }
#' 
#' \item{\code{src/ManifoldOptim/ManifoldOptim.cpp}:}{
#' Contains the main function \code{ManifoldOptim} which takes a problem
#' constructed in R, sets it up in ROPTLIB, runs it, and returns the result.
#' }
#' 
#' \item{\code{src/ManifoldOptim/ManifoldOptimModule.cpp}:}{
#' Defines an Rcpp module for ManifoldOptim which exposes C++ classes such as
#' \code{RProblem}. This module provides the most common means in which R users
#' will interact with ManifoldOptim.
#' }
#' 
#' \item{\code{src/ManifoldOptim/ManifoldFactory.h}:}{
#' The \code{GetManifold} function constructs a Manifold object based on its
#' name and dimensions. Manifold classes are defined in ROPTLIB.
#' }
#' 
#' \item{\code{src/ManifoldOptim/ProblemAdapter.h}:}{
#' Defines the \code{ProblemAdapter} class, which takes a
#' \code{ManifoldOptimProblem}, which is defined in the ManifoldOptim
#' API, and plugs it into the ROPTLIB API as an ROPTLIB \code{Problem}
#' subclass.
#' }
#' 
#' \item{\code{src/ManifoldOptim/RProblem.h}:}{
#' Defines the \code{RProblem} class, which allows the objective, gradient,
#' and Hessian functions to be defined in R. When a function in the ROPTLIB
#' library invokes the objective, gradient, or Hessian, this class invokes
#' the appropriate function in R.
#' }
#' 
#' \item{\code{src/ManifoldOptim/SolverFactory.h}:}{
#' The \code{GetSolver} function constructs a Solver object based on its name,
#' a given \code{Problem}, an initial value, and an initial Hessian. Solver
#' classes are defined in ROPTLIB.
#' }
#' 
#' \item{\code{src/ManifoldOptim/Util.h}:}{
#' Defines a few utility functions, especially to assist in translating between
#' the ManifoldOptim C++ API and the ROPTLIB API.
#' }
#' 
#' \item{\code{src/ManifoldOptim/VariableFactory.h}:}{
#' The \code{GetVariable} function returns an optimization variable suitable
#' for a given Manifold, based on its name and dimension. Optimization
#' variables for supported Manifolds are defined in ROPTLIB.
#' }
#'
#' \item{\code{inst/include/ManifoldOptimException.h}:}{
#' Defines \code{ManifoldOptimException}, which is a subclass of STL
#' \code{exception}.
#' }
#' 
#' \item{\code{inst/include/ManifoldOptim.h}:}{
#' For users of the ManifoldOptim C++ API, this is the main header file to
#' include. For an example, see \code{inst/examples/brockett/cpp_sourceCpp/}.
#' }
#' 
#' \item{\code{inst/include/ManifoldOptimProblem.h}:}{
#' Defines \code{ManifoldOptimProblem}, which is the base class
#' for all optimization problems in the ManifoldOptim API. This class
#' facilitates writing problems with Armadillo, which can be
#' instantiated and manipulated in R, and solved through ROPTLIB.
#' This class assumes only that the optimization variable is a
#' one-dimensional vector; the user must reshape it into the appropriate
#' form (e.g. a matrix or list of matrices) when evaluating the objective,
#' gradient, and Hessian functions.
#' }
#' }
#' @name Design of C++ code
#' @references 
#' Wen Huang, P.A. Absil, K.A. Gallivan, Paul Hand (2016a). "ROPTLIB: an
#' object-oriented C++ library for optimization on Riemannian manifolds."
#' Technical Report FSU16-14, Florida State University.
#' 
#' Conrad Sanderson and Ryan Curtin. Armadillo: a template-based C++ library
#' for linear algebra. Journal of Open Source Software, Vol. 1, pp. 26, 2016.
#'
#' S. Martin, A. Raim, W. Huang, and K. Adragni (2020). "ManifoldOptim: 
#' An R Interface to the ROPTLIB Library for Riemannian Manifold Optimization."
#' Journal of Statistical Software, 93(1):1-32.
NULL

