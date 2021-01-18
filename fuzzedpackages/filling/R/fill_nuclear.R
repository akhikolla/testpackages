#' Low-Rank Completion with Nuclear Norm Optimization
#'
#' In many circumstances, it is natural to assume that there exists an underlying
#' low-rank structure. The assumption of \emph{low-rank} property leads to an optimization problem
#' for matrix completion problem,
#' \deqn{\mathrm{minimize}\quad rank(X)}
#' \deqn{\mathrm{s.t}~~ X_{ij}=A_{ij} ~~\mathrm{for}~~ A_{ij} \in E}
#' where \eqn{A_{ij}\in E} means the \eqn{(i,j)}-th entry of data matrix \eqn{A} is not missing. The objective
#' function can be further relaxed by nuclear norm
#' \deqn{\|X\|_* = \sum \sigma_i(X)}
#' where \eqn{\sigma_i (X)} is \eqn{i}-th singular value of the matrix \eqn{X}. Note that
#' for modeling purpose, we adopted closeness parameter \code{tolerance} for equality constraint.
#' \pkg{CVXR} package was used in implementation. Computational efficiency may not be guaranteed for large data matrix.
#'
#'
#' @examples
#' \dontrun{
#' ## load image data of 'lena64'
#' data(lena64)
#'
#' ## transform 5% of entries into missing
#' A <- aux.rndmissing(lena64, x=0.05)
#'
#' ## apply the method
#' filled <- fill.nuclear(A)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' image(A, col=gray((0:100)/100), axes=FALSE, main="5% missing")
#' image(filled$X, col=gray((0:100)/100), axes=FALSE, main="processed")
#' par(opar)
#' }
#'
#' @param A an \eqn{(n\times p)} partially observed matrix.
#' @param tolerance level of tolerance for entrywise equality condition.
#'
#' @return a named list containing \describe{
#' \item{X}{an \eqn{(n\times p)} matrix after completion.}
#' \item{norm}{solution of the minimization problem; approximate rank.}
#' \item{cvxr.status}{``optimal'' denotes the problem was solved. See \code{\link[CVXR]{psolve}} for more details on solvability.}
#' \item{cvxr.niters}{the number of iterations taken.}
#' \item{cvxr.solver}{type of solver used by \pkg{CVXR}.}
#' }
#'
#' @references
#' \insertRef{candes_exact_2009}{filling}
#'
#' @rdname fill_nuclear
#' @export
fill.nuclear <- function(A, tolerance=1e-3){
  #-----------------------------------------------------------------
  ## PREPROCESSING
  #   1. data check and dimension
  M = check_data(A)
  n = nrow(M)
  p = ncol(M)
  #   2. index for non-missing entries
  NonMissing <- (!is.na(M))
  M[!NonMissing] = 0
  #   3. tolerance condition
  if ((check_na(tolerance))||(is.infinite(tolerance))||(tolerance < sqrt(.Machine$double.eps))){
    tolerance = sqrt(.Machine$double.eps)
  }
  #   4. full column or row
  if (check_bycol(M)==FALSE){   message("* fill.nuclear : there exists at least one column full of missing entries.")}
  if (check_bycol(t(M))==FALSE){message("* fill.nuclear : there exists at least one row full of missing entries.")}


  #-----------------------------------------------------------------
  ## CVXR Approach
  X   <- CVXR::Variable(n,p)
  obj <- CVXR::norm_nuc(X)

  masked_X = CVXR::multiply(NonMissing, X)
  masked_M = CVXR::multiply(NonMissing, M)

  # if (utils::packageVersion("CVXR") < "1.0.0"){
  #   masked_X = mul_elemwise(NonMissing, X)
  #   masked_M = mul_elemwise(NonMissing, M)
  # } else {
  #   masked_X = CVXR::multiply(NonMissing, X)
  #   masked_M = CVXR::multiply(NonMissing, M)
  # }

  abs_diff = abs(masked_X-masked_M)
  constr   = list(abs_diff <= tolerance)

  cvxrProb <- Problem(Minimize(obj),constr)
  cvxrresult <- solve(cvxrProb)


  #-----------------------------------------------------------------
  ## RETURN RESULTS
  result = list()
  result$X = cvxrresult$getValue(X)
  result$norm = cvxrresult$value
  result$cvxr.status = cvxrresult$status
  result$cvxr.niters = cvxrresult$num_iters
  result$cvxr.solver = cvxrresult$solver
  return(result)
}
