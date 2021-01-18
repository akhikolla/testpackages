# AUXILIARY FUNCTIONS -----------------------------------------------------
# 01.
# 02. aux.rndmissing


# 02. aux.rndmissing ------------------------------------------------------
#' Randomly assign NAs to the data matrix with probability \code{x}
#'
#' \code{aux.rndmissing} randomly selects \eqn{100\cdot x}\% of entries from
#' a given data matrix and turns them into missing entries, i.e., their values
#' become \code{NA}.
#'
#' @param A an \eqn{(n\times p)} data matrix.
#' @param x percentage of turning current entries into missing (\code{NA}).
#'
#' @return an \eqn{(n\times p)} data matrix with missing entries at proportion \eqn{x}.
#'
#'
#' @examples
#' # load lena64 image matrix
#' data(lena64)
#'
#' # generate 10% of missing values
#' lena64_miss <- aux.rndmissing(lena64)
#'
#' # visualize
#' par(mfrow=c(1,2))
#' image(lena64, axes=FALSE, main="original image")
#' image(lena64_miss, axes=FALSE, main="10% missing entries")
#'
#' @rdname aux_rndmissing
#' @export
aux.rndmissing <- function(A, x=0.10){
  # checking and coercing to be matrix.
  B = check_data(A)
  # check x
  x = as.double(x)
  if ((length(as.vector(x))!=1)||(x<=0)||(x>=1)){
    stop("* aux.rndmissing : percentage x should be in (0,1).")
  }
  # generating missing mechanism.
  nelem = nrow(B)*ncol(B)
  nmiss = ceiling(as.double(nelem)*x)
  vecid = c(rep(TRUE,nmiss), rep(FALSE,nelem-nmiss))
  vecid = sample(vecid, length(vecid), replace=TRUE)
  B[vecid] = NA
  # return output
  return(B)
}
