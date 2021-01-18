#' Format genotype matrix
#'
#' Convert a 4-column genotype matrix or data.frame into a simplified numerical format.
#'
#' @param X Genotype matrix or data.frame.
#' @param Y.only Logical.  Whether or not to output only \code{Y} (see return).
#' @return A list with elements
#' \itemize{
#' \code{Y}: A 4-column numeric matrix.  Each row is sorted with zeros indicating missing alleles.
#' \code{A}: Vector of allele names.
#' \code{G}: Unique allele combinations, i.e. \code{Y} without duplicated rows.
#' }
#' @keywords internal
geno.format <- function(X, Y.only = TRUE) {
  # convert to numeric
  Y <- as.matrix(X)
  if(!is.numeric(Y)) {
    if(any(is.na(Y)) && any(Y[!is.na(Y)] == "")) {
      stop("Cannot supply both empty character and NA.")
    }
    Y[Y == ""] <- NA
    A <- factor(Y)
    Y <- matrix(as.numeric(A), ncol = 4)
    A <- levels(A)
  } else {
    A <- levels(factor(Y))
  }
  Y[is.na(Y)] <- 0
  # sort numeric
  Y <- t(apply(Y, 1, function(x) {
    x <- sort(x)
    c(x, rep(0, 4-length(x)))
  }))
  # unique allele combinations observed
  G <- unique(colSort(Y))
  if(Y.only) {
    ans <- Y
  } else {
    ans <- list(Y = Y, A = A, G = G)
  }
  ans
}
