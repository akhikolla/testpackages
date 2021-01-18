#' HardImpute : Generalized Spectral Regularization
#'
#' If the assumed underlying model has sufficiently many zeros, the LASSO type
#' shrinkage estimator is known to overestimate the number of non-zero coefficients.
#' \code{fill.HardImpute} aims at overcoming such difficulty via low-rank assumption
#' and hard thresholding idea, well-known concept in conventional regression analysis.
#' In algorithmic aspect, it takes output of \code{SoftImpute} as warm-start matrices
#' for iterative estimation process.
#'
#' @param A an \eqn{(n\times p)} partially observed matrix.
#' @param lambdas a length-\eqn{t} vector regularization parameters.
#' @param maxiter maximum number of iterations to be performed.
#' @param tol stopping criterion for an incremental progress.
#' @param rk assumed rank of the matrix.
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
#' ## apply the method with 3 rank conditions
#' fill1 <- fill.HardImpute(A, lambdas=c(500,100,50), rk=10)
#' fill2 <- fill.HardImpute(A, lambdas=c(500,100,50), rk=50)
#' fill3 <- fill.HardImpute(A, lambdas=c(500,100,50), rk=100)
#'
#' ## visualize only the last ones from each run
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,2), pty="s")
#' image(A, col=gray((0:100)/100), axes=FALSE, main="5% missing")
#' image(fill1$X[,,3], col=gray((0:100)/100), axes=FALSE, main="Rank 10")
#' image(fill2$X[,,3], col=gray((0:100)/100), axes=FALSE, main="Rank 50")
#' image(fill3$X[,,3], col=gray((0:100)/100), axes=FALSE, main="Rank 100")
#' par(opar)
#'
#' @references
#' \insertRef{mazumder_spectral_2010}{filling}
#'
#' @seealso \code{\link{fill.SoftImpute}}
#' @rdname fill_HardImpute
#' @export
fill.HardImpute <- function(A, lambdas=c(10,1,0.1), maxiter=100, tol=1e-3, rk=(min(dim(A))-1)){
  #-----------------------------------------------------------------
  ## PREPROCESSING
  #   1. data check, dimension, and missing-1 matrix
  X = check_data(A)
  n = nrow(X)
  p = ncol(X)
  IDmat = (is.na(X))*1.0
  if (check_bycol(X)==FALSE){   message("* fill.HardImpute : there exists at least one column full of missing entries.")}
  if (check_bycol(t(X))==FALSE){message("* fill.HardImpute : there exists at least one row full of missing entries.")}
  #   2. lambda
  lambdas = sort(as.double(lambdas),decreasing=TRUE)
  if ((any(lambdas<=0))||(any(is.na(lambdas)))||(any(is.infinite(lambdas)))||(!is.vector(lambdas))){
    stop("* fill.HardImpute : 'lambdas' should be a vector of regularization parameters.")
  }
  #   3. maxiter
  maxiter = as.integer(maxiter)
  if ((length(as.vector(maxiter))!=1)||(maxiter<=1)||(is.na(maxiter))||(is.infinite(maxiter))){
    stop("* fill.HardImpute : 'lambda' should be a positive integer.")
  }
  #   4. tol
  tol = as.double(tol)
  if ((length(as.vector(tol))!=1)||(tol<=0)||(is.na(tol))||(is.infinite(tol))){
    stop("* fill.HardImpute : 'tol' should be a positive real number.")
  }
  #   5. rk : rank
  rk = as.integer(rk)
  if ((length(as.vector(rk))!=1)||(rk<1)||(rk>=min(n,p))){
    stop("* fill.HardImpute : 'k' should be a positive integer in [1,min(#columns,#rows)).")
  }

  #-----------------------------------------------------------------
  ## MAIN COMPUTATION
  #   1. run SoftImpute first.
  Minit = X
  Minit[is.na(Minit)] = mean(Minit[!is.na(Minit)])
  Msol = cpp_SoftImpute(X,IDmat,lambdas,tol,maxiter,Minit)
  #   2. run HardImpute
  Hsol = cpp_HardImpute(X,IDmat,lambdas,tol,maxiter,Msol,rk);

  #-----------------------------------------------------------------
  ## RETURN
  result = list()
  result$X = Hsol
  return(result)
}
