#' @useDynLib HypergeoMat
#' @importFrom Rcpp evalCpp
NULL

isPositiveInteger <- function(m){
  is.vector(m) && is.numeric(m) && length(m) == 1L && m > 0 && floor(m) == m
}

isSymmetricPositive <- function(M){
  isSymmetric(M) && all(eigen(M, symmetric = TRUE, only.values = TRUE)$values >= 0)
}

isNotNegativeInteger <- function(z){
  Im(z) != 0 || Re(z)>0 || Re(z) != trunc(Re(z))
}

isNumericOrComplex <- function(x){
  is.vector(x) && is.atomic(x) && (is.numeric(x) || is.complex(x))
}
