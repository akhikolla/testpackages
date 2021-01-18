#' Incomplete Gamma function of a matrix argument
#'
#' @description Evaluates the incomplete Gamma function of a matrix argument.
#'
#' @param m truncation weight of the summation, a positive integer
#' @param a real or complex parameter with \code{Re(a)>(p-1)/2}, where
#' \code{p} is the dimension (the order of the matrix)
#' @param x either a real or complex square matrix,
#' or a numeric or complex vector, the eigenvalues of the matrix
#'
#' @return A real or complex number.
#' @export
#'
#' @note This function is usually defined
#' for a symmetric real matrix or a Hermitian complex matrix.
#'
#' @references A. K. Gupta and D. K. Nagar.
#' \emph{Matrix variate distributions}. Chapman and Hall, 1999.
#'
#' @examples # for a scalar x, this is the incomplete Gamma function:
#' a <- 2
#' x <- 1.5
#' IncGamma(m = 15, a, x)
#' gsl::gamma_inc_P(a, x)
#' pgamma(x, shape = a, rate = 1)
IncGamma <- function(m, a, x){
  if(is.matrix(x)){
#    stopifnot(isSymmetric(x))
    p <- nrow(x)
  }else if(isNumericOrComplex(x)){
    stopifnot((p <- length(x)) > 0L)
  }else{
    stop("Invalid `x` argument")
  }
  stopifnot(
    isPositiveInteger(m),
    isNumericOrComplex(a),
    length(a) == 1L,
    Re(a) > (p-1)/2
  )
  if(is.matrix(x)){
    DET <- det(x)
  }else{
    DET <- prod(x)
  }
  DET^a * .mvbeta(a, (p+1)/2, p) * .hypergeomPFQ(m, a, a+(p+1)/2, -x, alpha = 2)
}
