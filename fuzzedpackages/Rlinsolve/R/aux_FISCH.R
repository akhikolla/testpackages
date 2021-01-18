#' Generate a 2-dimensional discrete Poisson matrix
#'
#' Poisson equation is one of most well-known elliptic partial differential equations. In order to
#' give a concrete example, a discrete Poisson matrix is generated, assuming we have \code{N} number of
#' grid points for each dimension under square domain. \emph{fisch} is a German word for Poisson.
#'
#' @param N the number of grid points for each direction.
#' @param sparse a logical; \code{TRUE} for returning sparse matrix, \code{FALSE} otherwise.
#'
#' @return an \eqn{(N^2 \times N^2)} matrix having block banded structure.
#' @examples
#' ## generate dense and sparse Poisson matrix of size 25 by 25.
#' A = aux.fisch(5, sparse=FALSE)
#' B = aux.fisch(5, sparse=TRUE)
#' (all(A==B)) # TRUE if two matrices are equal.
#'
#'
#' @references Golub, G. H. and Van Loan, C. F. (1996) \emph{Matrix Computations, 3rd Ed.}, pages 177–180.
#' @rdname aux_FISCH
#' @export
aux.fisch <- function(N,sparse=FALSE){
  if ((!is.numeric(N))||(is.na(N))||(is.infinite(N))||(length(N)>1)||(N<=2)){
    stop("* aux.fisch : an input 'N' should be a positive integer larger than 2.")
  }
  if (!is.logical(sparse)){
    stop("* aux.fisch : 'sparse' is a logical flag.")
  }
  siz  = as.integer(N)
  leng = siz*siz;
  dia  = matrix(0,nrow=siz,ncol=siz)

  diag(dia) = 4;
  for (i in 1:(siz-1)){
    dia[i,i+1] = -1;
    dia[i+1,i] = -1;
  }

  if (sparse){
    mat = Matrix(0,nrow=leng,ncol=leng)
  } else {
    mat = matrix(0,nrow=leng,ncol=leng)
  }
  if (!sparse){
    for (ib in 1:siz){
      mat[(1+(ib-1)*siz):(ib*siz),(1+(ib-1)*siz):(ib*siz)] = dia
    }
    for (ib in 1:(siz-1)){
      mat[(1+(ib-1)*siz):(ib*siz),(1+ib*siz):((ib+1)*siz)] = -diag(siz)
      mat[(1+ib*siz):((ib+1)*siz),(1+(ib-1)*siz):(ib*siz)] = -diag(siz)
    }
  } else {
    dia = Matrix(dia,sparse=TRUE)
    off = Matrix(-diag(siz),sparse=TRUE)
    for (ib in 1:siz){
      mat[(1+(ib-1)*siz):(ib*siz),(1+(ib-1)*siz):(ib*siz)] = dia
    }
    for (ib in 1:(siz-1)){
      mat[(1+(ib-1)*siz):(ib*siz),(1+ib*siz):((ib+1)*siz)] = off
      mat[(1+ib*siz):((ib+1)*siz),(1+(ib-1)*siz):(ib*siz)] = off
    }
  }
  return(mat)
}
