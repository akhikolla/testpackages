#' Distribution of the largest root
#'
#' Computes the cumulative distribution function of the largest root in the single and double
#' Wishart setting.
#'
#' If \eqn{S} follows a Wishart(\eqn{p,n}) distribution, e.g. if we can write \deqn{S = X^TX,} where
#' \eqn{X} is an \eqn{n x p} matrix with i.i.d rows coming from a \eqn{p}-variate standard normal,
#' then \code{singleWishart} gives the distribution of the largest root of \eqn{S}.
#'
#' As its name indicates, the double Wishart setting involves two Wishart variables: let \eqn{A} and
#' \eqn{B} be Wishart(\eqn{p,m}) and Wishart(\eqn{p,n}), respectively. If \eqn{A+B} is invertible,
#' then \code{doubleWishart} gives the distribution of the largest root of \deqn{(A+B)^-1B.}
#' Alternatively, it gives the distribution of the largest root of the determinental equation
#' \deqn{det(B - \theta(A+B)).}
#'
#' @param x Vector of numeric values at which to compute the CDF.
#' @param p,n,m Parameters of the single and double Wishart settings. See details.
#' @param type Character string. Select \code{type = "multi"} for multiprecision; select \code{type
#'   = "double"} for double precision. Defaults to adaptive selection of the precision type based
#'   on the input parameters.
#' @return Returns the value of the CDF at \code{x}.
#' @examples
#' x1 <- seq(0, 30, length.out = 50)
#' y1 <- singleWishart(x1, p = 5, n = 10)
#' plot(x1, y1, type='l')
#'
#' x2 <- seq(0, 1, length.out = 50)
#' y2 <- doubleWishart(x2, p = 10, n = 10, m = 200)
#' plot(x2, y2, type='l')
#' @export
#' @aliases doubleWishart singleWishart
#' @rdname largestRoot
singleWishart <- function(x, p, n, type = c("double", "multiple")) {
    # Check input
    stopifnot(all(x >= 0), p > 0, n > 0, isWhole(p), isWhole(n))

    n_min <- min(p, n)
    n_max <- max(p, n)

    if (missing(type)) {
        mprec <- precSingleWishart_bool(n_min, n_max)
        if (mprec) {
            message("Using multiprecision")
        } else {
            message("Using double precision")
        }
    } else {
        mprec <- switch(type,
                        "multiple" = TRUE,
                        "double" = FALSE,
                        NULL)
    }
    if (is.null(mprec)) {
        stop("Invalid value for argument \"type\"")
    }

    if (mprec) {
        result <- singleWishart_raw(x, n_min, n_max, mprec)
    } else {
        # There are some special cases
        specialCases <- list(c(2,2),
                             c(2,5),
                             c(3,3),
                             c(4,4))
        if (any(vapply(specialCases,
                      function(pair) all(pair == c(n_min, n_max)),
                      logical(1)))) {
            if (all(c(n_min, n_max) == c(2,2))) result <- F22(x)
            if (all(c(n_min, n_max) == c(2,5))) result <- F25(x)
            if (all(c(n_min, n_max) == c(3,3))) result <- F33(x)
            if (all(c(n_min, n_max) == c(4,4))) result <- F44(x)
        } else {
            result <- singleWishart_raw(x, n_min, n_max, FALSE)
        }
    }

    return(result)
}

#' @export
#' @rdname largestRoot
doubleWishart <- function(x, p, n, m, type = c("double", "multiple")) {
    # Check input
    stopifnot(all(x >= 0), all(x <= 1),
              p > 0, isWhole(p),
              n > 0, isWhole(n),
              m > 0, isWhole(m))

    if (missing(type)) {
        mprec <- precDoubleWishart_bool(p, n, m)
        if (mprec) {
            message("Using multiprecision")
        } else {
            message("Using double precision")
        }
    } else {
        mprec <- switch(type,
                        "multiple" = TRUE,
                        "double" = FALSE,
                        NULL)
    }
    if (is.null(mprec)) {
        stop("Invalid value for argument \"type\"")
    }

    # Convert to Chiani's notation
    sC <- p
    mC <- 0.5*(abs(n - p) - 1)
    nC <- 0.5*(abs(m - p) - 1)

    doubleWishart_raw(x, sC, mC, nC, mprec)
}

#' @useDynLib rootWishart, .registration=TRUE
#' @importFrom Rcpp sourceCpp
NULL
