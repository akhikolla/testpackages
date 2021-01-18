#' Sparse matrix
#' 
#' @description Constructs a sparse matrix, real or complex.
#'
#' @param i,j indices of the non-zero coefficients
#' @param Mij values of the non-zero coefficients; must be a vector of the same 
#'   length as \code{i} and \code{j} or a single number which will be recycled
#' @param nrows,ncols dimensions of the matrix
#' @param x a \code{SparseMatrix} object
#' @param ... ignored
#' @param M a matrix, real or complex
#'
#' @return A list with the class \code{SparseMatrix}.
#' @export
#' 
#' @name SparseMatrix
#'
#' @examples set.seed(666)
#' ( M <- matrix(rpois(50L, 1), 10L, 5L) )
#' asSparseMatrix(M)
SparseMatrix <- function(i, j, Mij, nrows, ncols){
  stopifnot(is.atomic(i), is.atomic(j), is.atomic(Mij))
  stopifnot(isStrictPositiveInteger(i), isStrictPositiveInteger(j))
  stopifnot(is.numeric(Mij) || is.complex(Mij))
  stopifnot(length(nrows) == 1L, length(ncols) == 1L)
  stopifnot(isStrictPositiveInteger(nrows), isStrictPositiveInteger(ncols))
  stopifnot(max(i) <= nrows, max(j) <= ncols)
  stopifnot(length(i) == length(j))
  stopifnot(length(Mij) == 1L || length(Mij) == length(i))
  if(length(Mij) == 1L){
    Mij <- rep(Mij, times = length(i))
  }
  out <- list(i = i-1L, j = j-1L, Mij = Mij, nrows = nrows, ncols = ncols)
  class(out) <- "SparseMatrix"
  out
}

#' @rdname SparseMatrix
#' @export
print.SparseMatrix <- function(x, ...){
  M <- matrix(".", nrow = x$nrows, ncol = x$ncols)
  Mij <- formatC(x$Mij)
  M[cbind(x$i+1L, x$j+1L)] <- Mij
  print(M, quote = FALSE)
}

#' @rdname SparseMatrix
#' @export
asSparseMatrix <- function(M){
  ij <- which(M != 0, arr.ind = TRUE)
  SparseMatrix(
    i = ij[,1L], j = ij[,2L], Mij = M[ij], nrows = nrow(M), ncols = ncol(M)
  )
}
