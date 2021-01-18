#' Imputation using Weighted K-nearest Neighbors
#'
#' One of the simplest idea to \emph{guess} missing entry is to use
#' portion of the data that has most similar characteristics across
#' all covariates. \code{fill.KNNimpute} follows such reasoning in that
#' it finds \eqn{K}-nearest neighbors based on observed variables and
#' uses weighted average of nearest elements to fill in the missing entry.
#' Note that when there are many missing entries, it's possible that there are
#' no \emph{surrogates} to be computed upon. Therefore, if there exists an entire
#' row or column full of missing entries, the algorithm stops.
#'
#' @param A an \eqn{(n\times p)} partially observed matrix.
#' @param k the number of neighbors to use.
#'
#' @return a named list containing \describe{
#' \item{X}{an \eqn{(n\times p)} matrix after completion.}
#' }
#'
#' @examples
#' ## load image data of 'lena128'
#' data(lena128)
#'
#' ## transform 5% of entries into missing
#' set.seed(5)
#' A <- aux.rndmissing(lena128, x=0.05)
#'
#' ## apply the method with 3 different neighborhood size
#' fill1 <- fill.KNNimpute(A, k=5)
#' fill2 <- fill.KNNimpute(A, k=25)
#' fill3 <- fill.KNNimpute(A, k=50)
#'
#' ## visualize only the last ones from each run
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,2), pty="s")
#' image(A, col=gray((0:100)/100), axes=FALSE, main="5% missing")
#' image(fill1$X, col=gray((0:100)/100), axes=FALSE, main="5-neighbor")
#' image(fill2$X, col=gray((0:100)/100), axes=FALSE, main="25-neighbor")
#' image(fill3$X, col=gray((0:100)/100), axes=FALSE, main="50-neighbor")
#' par(opar)
#'
#' @references
#' \insertRef{troyanskaya_missing_2001}{filling}
#'
#' @seealso \code{\link{fill.SVDimpute}}
#' @rdname fill_KNNimpute
#' @export
fill.KNNimpute <- function(A, k=ceiling(nrow(A)/2)){
  #-----------------------------------------------------------------
  ## PREPROCESSING
  #   1. data check, dimension, and missing-1 matrix
  X = check_data(A)
  N = nrow(X)
  P = ncol(X)
  if (check_bycol(X)==FALSE){   stop("* fill.KNNimpute : there exists at least one column full of missing entries.")}
  if (check_bycol(t(X))==FALSE){stop("* fill.KNNimpute : there exists at least one row full of missing entries.")}
  #   2. k : size of neighbor
  k = as.integer(k)
  if ((length(as.vector(k))!=1)||(k<1)||(k>=P)){
    stop("* fill.KNNimpute : 'k' should be a positive integer in [1,#(number of columns)).")
  }

  #-----------------------------------------------------------------
  ## MAIN COMPUTATION
  #   1. find missing-entry indices
  idxmiss = which(is.na(X), arr.ind=TRUE)
  if (length(nrow(idxmiss))==0){stop("* fill.KNNimpute : there is no missing entries.")}

  #   2. random permutation of rows
  rowperm = sample(1:nrow(idxmiss), replace=FALSE)
  idxmiss = idxmiss[rowperm,]

  #   3. iterate!
  for (i in 1:nrow(idxmiss)){
    # 3-1. current missing index
    crow = idxmiss[i,1]
    ccol = idxmiss[i,2]
    # 3-2. non-missings in
    obsrows = which(!is.na(X[,ccol])) # this is target for knn
    # 3-3. let's compute it.
    X[crow,ccol] = KNNimpute_single(X[crow,], X[obsrows, ], ccol, min(k,length(obsrows)))
  }

  #-----------------------------------------------------------------
  ## RETURN
  result = list()
  result$X = X
  return(result)
}


#' @keywords internal
#' @noRd
KNNimpute_single <- function(rowvec, mat, colidx, park){
  col_matmiss = which(is.na(colSums(mat)))
  col_vecmiss = which(is.na(rowvec))
  all_missing = union(col_matmiss, col_vecmiss)

  if (length(all_missing)==ncol(mat)){
    solution = mean(mat[,colidx])
  } else {
    newvec = rowvec[-all_missing]
    newmat = mat[,-all_missing]

    newknn = nabor::knn(newmat, query=matrix(newvec,ncol=length(newvec)), park)

    knnidx    = newknn$nn.idx
    knndistinv = 1/newknn$nn.dists
    knnweight = knndistinv/sum(knndistinv)

    solution = sum(as.vector(mat[knnidx, colidx])*knnweight)
  }
  return(solution)
}
