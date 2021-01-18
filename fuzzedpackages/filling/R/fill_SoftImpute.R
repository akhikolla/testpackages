#' SoftImpute : Spectral Regularization
#'
#' \code{fill.SoftImpute} implements convex relaxation techniques to generate
#' a sequence of regularized low-rank solutions for matrix completion problems.
#' For the nuclear norm optimization problem, it uses soft thresholding technique
#' iteratively in that the algorithm returns several matrices in accordance with
#' the provided vector of regularization parameters \eqn{\lambda} (\code{lambdas}).
#'
#' @param A an \eqn{(n\times p)} partially observed matrix.
#' @param lambdas a length-\eqn{t} vector regularization parameters.
#' @param maxiter maximum number of iterations to be performed.
#' @param tol stopping criterion for an incremental progress.
#'
#' @return a named list containing \describe{
#' \item{X}{an \eqn{(n\times p\times t)} cubic array after completion at each \code{lambda} value.}
#' }
#'
#' @examples
#' ## load image data of 'lena128'
#' data(lena128)
#'
#' ## transform 5% of entries into missing
#' A <- aux.rndmissing(lena128, x=0.05)
#'
#' ## apply the method with 3 lambda values
#' fill <- fill.SoftImpute(A, lambdas=c(500,100,50))
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,2), pty="s")
#' image(A, col=gray((0:100)/100), axes=FALSE, main="5% missing")
#' image(fill$X[,,1], col=gray((0:100)/100), axes=FALSE, main="lambda=500")
#' image(fill$X[,,2], col=gray((0:100)/100), axes=FALSE, main="lambda=100")
#' image(fill$X[,,3], col=gray((0:100)/100), axes=FALSE, main="lambda=50")
#' par(opar)
#'
#' @references
#' \insertRef{mazumder_spectral_2010}{filling}
#'
#' @seealso \code{\link{fill.HardImpute}}
#' @rdname fill_SoftImpute
#' @export
fill.SoftImpute <- function(A, lambdas=c(10,1,0.1), maxiter=100, tol=1e-3){
  #-----------------------------------------------------------------
  ## PREPROCESSING
  #   1. data check, dimension, and missing-1 matrix
  X = check_data(A)
  n = nrow(X)
  p = ncol(X)
  IDmat = (is.na(X))*1.0
  if (check_bycol(X)==FALSE){   message("* fill.SoftImpute : there exists at least one column full of missing entries.")}
  if (check_bycol(t(X))==FALSE){message("* fill.SoftImpute : there exists at least one row full of missing entries.")}
  #   2. lambda
  lambdas = sort(as.double(lambdas),decreasing=TRUE)
  if ((any(lambdas<=0))||(any(is.na(lambdas)))||(any(is.infinite(lambdas)))||(!is.vector(lambdas))){
    stop("* fill.SoftImpute : 'lambdas' should be a vector of regularization parameters.")
  }
  #   3. maxiter
  maxiter = as.integer(maxiter)
  if ((length(as.vector(maxiter))!=1)||(maxiter<=1)||(is.na(maxiter))||(is.infinite(maxiter))){
    stop("* fill.SoftImpute : 'lambda' should be a positive integer.")
  }
  #   4. tol
  tol = as.double(tol)
  if ((length(as.vector(tol))!=1)||(tol<=0)||(is.na(tol))||(is.infinite(tol))){
    stop("* fill.SoftImpute : 'tol' should be a positive real number.")
  }

  #-----------------------------------------------------------------
  ## MAIN COMPUTATION
  Minit = X
  Minit[is.na(Minit)] = mean(Minit[!is.na(Minit)])
  Msol = cpp_SoftImpute(X,IDmat,lambdas,tol,maxiter,Minit)

  #-----------------------------------------------------------------
  ## RETURN
  result = list()
  result$X = Msol
  return(result)
}
