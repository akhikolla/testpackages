#' Iterative Regression against Right Singular Vectors
#'
#' Singular Value Decomposition (SVD) is the best low-rank approximation of a given matrix.
#' \code{fill.SVDimpute} exploits such idea. First, it starts with simple filling using
#' column mean values for filling. Second, it finds SVD of a current matrix. Then,
#' each row vector is regressed upon top-\eqn{k} right singular vectors. Missing entries are
#' then filled with predicted estimates.
#'
#' @param A an \eqn{(n\times p)} partially observed matrix.
#' @param k the number of regressors to be used.
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
#' set.seed(5)
#' A <- aux.rndmissing(lena128, x=0.05)
#'
#' ## apply the method with 3 different number of regressors
#' fill1 <- fill.SVDimpute(A, k=5)
#' fill2 <- fill.SVDimpute(A, k=25)
#' fill3 <- fill.SVDimpute(A, k=50)
#'
#' ## visualize only the last ones from each run
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,2), pty="s")
#' image(A, col=gray((0:100)/100), axes=FALSE, main="5% missing")
#' image(fill1$X, col=gray((0:100)/100), axes=FALSE, main="5 regressors")
#' image(fill2$X, col=gray((0:100)/100), axes=FALSE, main="25 regressors")
#' image(fill3$X, col=gray((0:100)/100), axes=FALSE, main="50 regressors")
#' par(opar)
#' }
#'
#'
#' @references
#' \insertRef{troyanskaya_missing_2001}{filling}
#'
#' @seealso \code{\link{fill.KNNimpute}}
#' @rdname fill_SVDimpute
#' @export
fill.SVDimpute <- function(A, k=ceiling(ncol(A)/2), maxiter=100, tol=1e-2){
  #-----------------------------------------------------------------
  ## PREPROCESSING
  #   1. data check, dimension, and missing-1 matrix
  X = check_data(A)
  N = nrow(X)
  P = ncol(X)
  if (check_bycol(X)==FALSE){message("* fill.SVDimpute : there exists at least one column full of missing entries.")}
  if (check_bycol(t(X))==FALSE){stop("* fill.SVDimpute : there exists at least one row full of missing entries.")}
  #   2. k : number of top singular vectors to be used
  k = as.integer(k)
  if ((length(as.vector(k))!=1)||(k<1)||(k>P)){
    stop("* fill.SVDimpute : 'k' is the number of right singular vectors to be used, [1.#(columns)).")
  }
  #   3. maxiter
  maxiter = as.integer(maxiter)
  if ((length(as.vector(maxiter))!=1)||(maxiter<=1)||(is.na(maxiter))||(is.infinite(maxiter))){
    stop("* fill.SVDImpute : 'lambda' should be a positive integer.")
  }
  #   4. tol
  tol = as.double(tol)
  if ((length(as.vector(tol))!=1)||(tol<=0)||(is.na(tol))||(is.infinite(tol))){
    stop("* fill.SVDImpute : 'tol' should be a positive real number.")
  }

  #-----------------------------------------------------------------
  ## MAIN COMPUTATION
  #   1. get index for each row of NAs.
  idxmiss = list()
  for (i in 1:N){
    idxmiss[[i]] = which(is.na(X[i,]))
  }
  #   2. fill with simple 'mean' method.
  Xold = fill.simple(X, method="mean")$X
  Xnew = array(0,c(N,P))
  #   3. Iterate
  inctol = 10000.0
  iter   = 0;
  while (inctol > tol){
    # 3-1. compute top-k SVD
    Vk = RSpectra::svds(Xold, k)$v
    # 3-2. fill in
    Xnew = SVDimpute_single(Xold, Vk, idxmiss, N, P)
    # 3-3. compute reltol
    inctol = norm(Xnew-Xold,"f")/norm(Xold,"f")
    # 3-4. update
    Xold = Xnew
    # 3-5. iteration
    iter = iter + 1
    if (iter >= maxiter){
      break
    }
  }
  #-----------------------------------------------------------------
  ## RETURN
  result = list()
  result$X = Xold
  return(result)
}


#' @keywords internal
#' @noRd
SVDimpute_single <- function(Mold, V, idx, n, p){
  Mnew = Mold
  for (i in 1:nrow(Mold)){
    vecLHS = Mold[i,-idx[[i]]]
    matRHS = V[-idx[[i]],]
    solvec = as.vector(lm(vecLHS~-1+matRHS)$coefficients)
    estLHS = as.vector(V%*%matrix(solvec, ncol=1))

    Mnew[i,idx[[i]]]  = estLHS[idx[[i]]] # missing
  }
  return(Mnew)
}
