#-----------------------------------------------------------------------#
# Package: PaRametric sImplex Method for spArse Learning                #
# ParametricSimplexMethod_solver() : Solve given dantzig selector problem    #
# in parametric simplex method                                          #
#-----------------------------------------------------------------------#

#' Solve given problem in parametric simplex method
#'
#' @param A \code{A} is an \code{n} by \code{d} data matrix
#' @param b \code{b} is a length \code{n} response vector
#' @param b_bar \code{b_bar} is a length \code{n} vector time to parameter in constraints.
#' @param c \code{c} is a length \code{d} vector in target function.
#' @param c_bar \code{c_bar} is a length \code{d} vector time to parameter in target function
#' @param B_init \code{B_init} is the index of initial basic colume.
#' @param max_it This is the number of the maximum path length one would like to achieve. The default length is \code{50}.
#' @param lambda_threshold The parametric simplex method will stop when the calculated parameter is smaller than lambda. The default value is \code{0.01}.
#' @return
#' An object with S3 class \code{"primal"} is returned:
#' \item{data}{
#'   The \code{n} by \code{d} data matrix from the input
#' }
#' \item{response}{
#'   The length \code{n} response vector from the input
#' }
#' \item{beta}{
#'   A matrix of regression estimates whose columns correspond to regularization parameters for parametric simplex method.
#' }
#' \item{beta0}{
#'   A vector of regression estimates whose index correspond to regularization parameters for parametric simplex method.
#' }
#' \item{df}{
#'   The degree of freecom (number of nonzero coefficients) along the solution path.
#' }
#' \item{value}{
#'   The sequence of optimal value of the object function corresponded to the sequence of lambda.
#' }
#' \item{iterN}{
#'   The number of iteration in the program.
#' }
#' \item{lambda}{
#'   The sequence of regularization parameters \code{lambda} obtained in the program.
#' }
#' \item{type}{
#'   The type of the problem, such as \code{Dantzig} and \code{SparseSVM}.
#' }
#' @examples
#' ## This example show how to use PSM_solver() to solve dantzig problem.
#' ## Generate the design matrix and coefficient vector
#' n = 100 # sample number
#' d = 250 # sample dimension
#' c = 0.5 # correlation parameter
#' s = 20  # support size of coefficient
#' set.seed(1024)
#' X = scale(matrix(rnorm(n*d),n,d)+c*rnorm(n))/sqrt(n-1)*sqrt(n)
#' beta = c(rnorm(s), rep(0, d-s))
#' ## Generate response using Gaussian noise, and solve the solution path
#' noise = rnorm(n)
#' Y = X%*%beta + noise
#' ## Define parameters for dantzig problem
#' XtX = t(X)%*%X
#' A = cbind(cbind(rbind(XtX,-XtX),-rbind(XtX,-XtX)),diag(rep(1,2*d)))
#' b = rbind(t(X)%*%Y,-t(X)%*%Y)
#' c = c(rep(-1,2*d),rep(0,2*d))
#' c_bar = rep(0,4*d)
#' b_bar = rep(1,2*d)
#' B_init = seq(2*d,4*d-1)
#' ## Dantzig selection solved with parametric simplex method
#' fit.dantzig = PSM_solver(A, b, b_bar, c, c_bar, B_init, max_it = 50, lambda_threshold = 0.01)
#' ###lambdas used
#' print(fit.dantzig$lambda)
#' ## number of nonzero coefficients for each lambda
#' print(fit.dantzig$df)
#' ## Visualize the solution path
#' plot(fit.dantzig)
#' @seealso \code{\link{primal-package}}
#' @export
PSM_solver <- function(A, b, b_bar, c, c_bar, B_init = NULL, max_it = 50, lambda_threshold = 0.01) {
  begt <- Sys.time()
  n <- nrow(A)
  d <- ncol(A)
  t <- 0
  lambdalist <- rep(0, max_it)
  x_list <- matrix(0, d, max_it)
  y_list <- rep(0, max_it)
  if(is.null(B_init)){
    message("The initial Basic is needed!")
  }
  str <- .C("ParametricSimplexMethod_api", as.integer(n), as.integer(d), as.double(t(A)), as.double(b), as.double(b_bar), as.double(c),
            as.double(c_bar),as.integer(max_it), as.double(lambda_threshold), as.integer(t),
            as.double(lambdalist), as.double(x_list), as.double(y_list),as.integer(B_init), PACKAGE = "PRIMAL")
  t <- unlist(str[10])
  x_list <- matrix(unlist(str[12])[1:(d * t)], d, t)
  df <- c()
  for (i in 1:t) {
    df[i] <- sum(x_list[, i] != 0)
  }
  runt <- Sys.time() - begt
  ans <- list(type = "PSM",
              data = A,
              response = b,
              beta = Matrix(x_list),
              df = df,
              value = unlist(str[13])[1:t],
              iterN = t,
              lambda = unlist(str[11])[1:t],
              runtime = runt)
  class(ans) <- "primal"
  return(ans)
}


