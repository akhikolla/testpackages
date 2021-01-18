#' OptSpace
#'
#' OptSpace is an algorithm for matrix completion when a matrix is partially observed. It
#' performs what authors called \emph{trimming} and \emph{projection} repeatedly based on
#' singular value decompositions. Original implementation is borrowed from \pkg{ROptSpace} package,
#' which was independently developed by the maintainer. See \code{\link[ROptSpace]{OptSpace}} for more details.
#'
#' @examples
#' \dontrun{
#' ## load image data of 'lena64'
#' data(lena64)
#'
#' ## transform 5% of entries into missing
#' A <- aux.rndmissing(lena64, x=0.05)
#'
#' ## apply the method with different rank assumptions
#' filled10 <- fill.OptSpace(A, ropt=10)
#' filled20 <- fill.OptSpace(A, ropt=20)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' image(A, col=gray((0:100)/100), axes=FALSE, main="5% missing")
#' image(filled10$X, col=gray((0:100)/100), axes=FALSE, main="rank 10")
#' image(filled20$X, col=gray((0:100)/100), axes=FALSE, main="rank 20")
#' par(opar)
#' }
#'
#'
#' @param A an \eqn{(n\times p)} partially observed matrix.
#' @param ropt \code{NA} to guess the rank, or a positive integer as a pre-defined rank.
#' @param niter maximum number of iterations allowed.
#' @param tol stopping criterion for reconstruction in Frobenius norm.
#'
#' @return a named list containing \describe{
#' \item{X}{an \eqn{(n\times p)} matrix after completion.}
#' \item{error}{a vector of reconstruction errors for each successive iteration.}
#' }
#'
#' @references
#' \insertRef{keshavan_matrix_2010}{filling}
#'
#' @rdname fill_OptSpace
#' @export
fill.OptSpace <- function(A, ropt=NA, niter=50, tol=1e-6){
  #-----------------------------------------------------------------
  ## PREPROCESSING
  #   1. data check and dimension
  M = check_data(A)
  n = nrow(M)
  p = ncol(M)
  #   2. full column or row
  if (check_bycol(M)==FALSE){   message("* fill.optspace : there exists at least one column full of missing entries.")}
  if (check_bycol(t(M))==FALSE){message("* fill.optspace : there exists at least one row full of missing entries.")}

  #-----------------------------------------------------------------
  ## ROptSpace
  fill_ropt = ropt
  fill_niter= niter
  fill_tol  = tol
  output = ROptSpace::OptSpace(M, ropt=fill_ropt, niter=fill_niter, tol=fill_tol, showprogress=FALSE)


  #-----------------------------------------------------------------
  ## RETURN RESULTS
  result = list()
  result$X = (output$X %*% output$S %*% t(output$Y))
  result$error = output$dist
  return(result)
}



