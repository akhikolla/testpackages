#' Singular Value Thresholding for Nuclear Norm Optimization
#'
#' \code{fill.SVT} is an iterative updating scheme for Nuclear Norm Minimization problem. An unconstrained
#' parahrase of the problem introduced in \code{\link{fill.nuclear}} is
#' \deqn{\mathrm{minimize}\quad \frac{1}{2}\|P_{\Omega}(X-A) \|_F^2 + \lambda \| X \|_*}
#' where \eqn{P_{\Omega}(X)=X_{ij}} if it is observed, or \eqn{0} otherwise.
#' It performs iterative shrinkage on newly computed singular values.
#'
#' @param A an \eqn{(n\times p)} partially observed matrix.
#' @param lambda a regularization parameter.
#' @param maxiter maximum number of iterations to be performed.
#' @param tol stopping criterion for an incremental progress.
#'
#' @return a named list containing \describe{
#' \item{X}{an \eqn{(n\times p)} matrix after completion.}
#' }
#'
#' @examples
#' \dontrun{
#' ## load image data of 'lena128'
#' data(lena128)
#'
#' ## transform 5% of entries into missing
#' A <- aux.rndmissing(lena128, x=0.05)
#'
#' ## apply the method
#' fill1 <- fill.SVT(A, lambda=0.1)
#' fill2 <- fill.SVT(A, lambda=1.0)
#' fill3 <- fill.SVT(A, lambda=20)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,2), pty="s")
#' image(A, col=gray((0:100)/100), axes=FALSE, main="5% missing")
#' image(fill1$X, col=gray((0:100)/100), axes=FALSE, main="lbd=0.1")
#' image(fill2$X, col=gray((0:100)/100), axes=FALSE, main="lbd=1")
#' image(fill3$X, col=gray((0:100)/100), axes=FALSE, main="lbd=10")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{cai_singular_2010}{filling}
#'
#' @seealso \code{\link{fill.nuclear}}
#' @rdname fill_SVT
#' @export
fill.SVT <- function(A, lambda=1.0, maxiter=100, tol=1e-3){
  #-----------------------------------------------------------------
  ## PREPROCESSING
  #   1. data check, dimension, and missing-1 matrix
  X = check_data(A)
  n = nrow(X)
  p = ncol(X)
  IDmat = (is.na(X))*1.0
  if (check_bycol(X)==FALSE){   message("* fill.SVT : there exists at least one column full of missing entries.")}
  if (check_bycol(t(X))==FALSE){message("* fill.SVT : there exists at least one row full of missing entries.")}
  #   2. lambda
  lambda = as.double(lambda)
  if ((length(as.vector(lambda))!=1)||(lambda<=0)||(is.na(lambda))||(is.infinite(lambda))){
    stop("* fill.SVT : 'lambda' should be a positive real number.")
  }
  #   3. maxiter
  maxiter = as.integer(maxiter)
  if ((length(as.vector(maxiter))!=1)||(maxiter<=1)||(is.na(maxiter))||(is.infinite(maxiter))){
    stop("* fill.SVT : 'lambda' should be a positive integer.")
  }
  #   4. tol
  tol = as.double(tol)
  if ((length(as.vector(tol))!=1)||(tol<=0)||(is.na(tol))||(is.infinite(tol))){
    stop("* fill.SVT : 'tol' should be a positive real number.")
  }

  #-----------------------------------------------------------------
  ## MAIN COMPUTATION
  #   1. initial values for M : fill with random numbers
  Minit = X
  Minit[is.na(Minit)] = mean(Minit[!is.na(Minit)])

  #   2. main compute !
  Msol = cpp_nSVD(X, IDmat, Minit, lambda, tol, maxiter)

  #-----------------------------------------------------------------
  ## RETURN
  result = list()
  result$X = Msol
  return(result)
}


