#' Check if a numeric value
#'
#' this function checks if the input is a scalar (double)
#' @param x is the input
#' @param n is the input length
#' @keywords internal

is_numeric <- function(x, n) {
    (typeof(x) == "double" || typeof(x) == "integer")  && length(x) == n && all(!is.na(x))
}

#' Check if a positive numeric value
#'
#'  this function checks if the input is a positive scalar (double)
#' @param x is the input
#' @param n is the input length
#' @keywords internal

is_positive_numeric <- function(x, n) {
    is_numeric(x, n) && all(x > 0)
}

#' Check if a numeric vector of length n
#'
#' this function checks if the input is a numeric (double) vector
#' @param x is the input
#' @param n is the input vector length
#' @keywords internal

is_numeric_vector <- function(x, n) {
    is_numeric(x, n) && is.vector(x)
}

#' Check if a numeric matrix of dimension \eqn{n \times m}{n x m}
#'
#' this function checks if the input is a numeric (double) matrix of dimension \eqn{n \times m}{n x m}
#' @param x is the input
#' @param n is the input matrix rows
#' @param m is the input matrix columns
#' @keywords internal

is_numeric_matrix <- function(x, n, m) {
    is_numeric(x, n * m) && is.matrix(x) && all(dim(x) == c(n, m))
}

#' Check if a symmetric positive definite numeric matrix of dimension \eqn{n \times n}{n x n}
#'
#' this function checks if the input is a symmetric positive definite matrix
#' @param x is the input
#' @param n is the input matrix dimension (assuming a square matrix)
#' @keywords internal

is_sympd_matrix <- function(x, n) {
    # if(!is.matrix(x))
    #     stop("must be a square matrix with number of rows = number of columns")
    # if (nrow(x) != ncol(x))
    #     stop("must be a square matrix with number of rows = number of columns")
    is_numeric_matrix(x, n, n) && n == n && isSymmetric(x) && all(eigen(x)$values > 0)
}

#' Check if value is an integer or integer-like
#'
#' this function checks if the input is an integer scalar (integer-like value -- i.e., both 1L and 1.0 pass this check)
#' @param x is the input
#' @param n is the number of inputs
#' @keywords internal

is_integer <- function(x, n) {
    if (is_numeric(x, n) && length(x) == n) {
        if (all(x == as.integer(x))) {
            return(TRUE)
        } else {
            return(FALSE)
        }
    } else {
        return(FALSE)
    }
    # typeof(x) == "integer" && length(x) == n
}

#'  Check if an integer matrix of dimension \eqn{n \times m}{n x m}
#'
#' this function checks if the input is an integer matrix of dimension \eqn{n \times m}{n x m}
#' @param x is the input
#' @param n is the input matrix rows
#' @param m is the input matrix columns
#' @keywords internal

is_integer_matrix <- function(x, n, m) {
    is_integer(x, n * m) && is.matrix(x) && all(dim(x) == c(n, m))
}



#' Check if value is a positive integer or integer-like
#'
#' this function checks if the input is a positive integer scalar (integer-like value -- i.e., both 1L and 1.0 pass this check)
#' @param x is the input
#' @param n is the number of inputs
#' @keywords internal

is_positive_integer <- function(x, n) {
    is_integer(x, n) && all(x > 0)
}






