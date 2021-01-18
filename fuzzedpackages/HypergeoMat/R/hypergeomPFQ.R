.hypergeomPFQ <- function(m, a, b, x, alpha){
  if(is.matrix(x)){
    x <- eigen(x, only.values = TRUE)$values
  }
  if(is.null(a)){
    a <- numeric(0L)
  }
  if(is.null(b)){
    b <- numeric(0L)
  }
  if(is.complex(a) || is.complex(b)){
    if(is.complex(x)){
      hypergeom_Cplx_Cplx(m = m, a = a, b = b, x = x, alpha = alpha)
    }else{
      hypergeom_Cplx_R(m = m, a = a, b = b, x = x, alpha = alpha)
    }
  }else{
    if(is.complex(x)){
      hypergeom_R_Cplx(m = m, a = a, b = b, x = x, alpha = alpha)
    }else{
      hypergeom_R_R(m = m, a = a, b = b, x = x, alpha = alpha)
    }
  }
}

#' Hypergeometric function of a matrix argument
#'
#' @description Evaluates a truncated hypergeometric function of a matrix
#' argument.
#'
#' @param m truncation weight of the summation, a positive integer
#' @param a the "upper" parameters, a numeric or complex vector,
#' possibly empty (or \code{NULL})
#' @param b the "lower" parameters, a numeric or complex vector,
#' possibly empty (or \code{NULL})
#' @param x either a real or complex square matrix,
#' or a numeric or complex vector, the eigenvalues of the matrix
#' @param alpha the alpha parameter, a positive number
#'
#' @return A real or a complex number.
#' @export
#'
#' @note The hypergeometric function of a matrix argument is usually defined
#' for a symmetric real matrix or a Hermitian complex matrix.
#'
#' @details This is an implementation of Koev & Edelman's algorithm
#' (see the reference). This algorithm is split into two parts: the case of
#' a scalar matrix (multiple of an identity matrix) and the general case.
#' The case of a scalar matrix is much faster (try e.g. \code{x = c(1,1,1)} vs
#' \code{x = c(1,1,0.999)}).
#'
#' @references Plamen Koev and Alan Edelman.
#' \emph{The Efficient Evaluation of the Hypergeometric Function of a Matrix Argument}.
#' Mathematics of Computation, 75, 833-846, 2006.
#'
#' @examples # a scalar x example, the Gauss hypergeometric function
#' hypergeomPFQ(m = 10, a = c(1,2), b = c(3), x = 0.2)
#' gsl::hyperg_2F1(1, 2, 3, 0.2)
#' # 0F0 is the exponential of the trace
#' X <- toeplitz(c(3,2,1))/10
#' hypergeomPFQ(m = 10, a = NULL, b = NULL, x = X)
#' exp(sum(diag(X)))
#' # 1F0 is det(I-X)^(-a)
#' X <- toeplitz(c(3,2,1))/100
#' hypergeomPFQ(m = 10, a = 3, b = NULL, x = X)
#' det(diag(3)-X)^(-3)
#' # Herz's relation for 1F1
#' hypergeomPFQ(m = 10, a = 2, b = 3, x = X)
#' exp(sum(diag(X))) * hypergeomPFQ(m = 10, a = 3-2, b = 3, x = -X)
#' # Herz's relation for 2F1
#' hypergeomPFQ(10, a = c(1,2), b = 3, x = X)
#' det(diag(3)-X)^(-2) *
#'   hypergeomPFQ(10, a = c(3-1,2), b = 3, -X %*% solve(diag(3)-X))
hypergeomPFQ <- function(m, a, b, x, alpha = 2){
  stopifnot(
    isPositiveInteger(m), 
    is.null(a) || isNumericOrComplex(a),
    is.null(b) || isNumericOrComplex(b),
    is.matrix(x) || isNumericOrComplex(x),
    length(x) > 0L, 
    is.vector(alpha) && is.atomic(alpha),
    length(alpha) == 1L,
    is.numeric(alpha),
    alpha > 0
  )
  .hypergeomPFQ(m = m, a = a, b = b, x = x, alpha = alpha)
}
