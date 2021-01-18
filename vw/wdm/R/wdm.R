#' Weighted Dependence Measures
#'
#' Computes a (possibly weighted) dependence measure between `x` and `y` if
#' these are vectors. If `x` and `y` are matrices then the measure between the
#' columns of `x` and the columns of `y` are computed.
#'
#' @param x a numeric vector, matrix or data frame.
#' @param y `NULL` (default) or a vector, matrix or data frame with compatible
#'   dimensions to x. The default is equivalent to `y = x`` (but more
#'   efficient).
#' @param method the dependence measure; see *Details* for possible values.
#' @param weights an optional vector of weights for the observations.
#' @param remove_missing if `TRUE`, all (pairswise) incomplete observations are
#'   removed; if `FALSE`, the function throws an error if there are incomplete
#'   observations.
#'
#' @details Available methods:
#' - `"pearson"`: Pearson correlation
#' - `"spearman"`: Spearman's \eqn{\rho}
#' - `"kendall"`: Kendall's \eqn{\tau}
#' - `"blomqvist"`: Blomqvist's \eqn{\beta}
#' - `"hoeffding"`: Hoeffding's \eqn{D}
#' Partial matching of method names is enabled.
#'
#' Spearman's \eqn{\rho} and Kendall's \eqn{\tau} are corrected for ties if
#' there are any.
#'
#' @export
#'
#' @examples
#' ##  dependence between two vectors
#' x <- rnorm(100)
#' y <- rpois(100, 1)  # all but Hoeffding's D can handle ties
#' w <- runif(100)
#' wdm(x, y, method = "kendall")               # unweighted
#' wdm(x, y, method = "kendall", weights = w)  # weighted
#'
#' ##  dependence in a matrix
#' x <- matrix(rnorm(100 * 3), 100, 3)
#' wdm(x, method = "spearman")               # unweighted
#' wdm(x, method = "spearman", weights = w)  # weighted
#'
#' ##  dependence between columns of two matrices
#' y <- matrix(rnorm(100 * 2), 100, 2)
#' wdm(x, y, method = "hoeffding")               # unweighted
#' wdm(x, y, method = "hoeffding", weights = w)  # weighted
#'
wdm <- function(x, y = NULL, method = "pearson", weights = NULL,
                remove_missing = TRUE) {
    ## preprocessing of arguments
    if (is.null(weights))
        weights <- numeric(0)
    if (is.data.frame(y))
        y <- as.matrix(y)
    if (is.data.frame(x))
        x <- as.matrix(x)
    check_wdm_inputs(x, y, weights, remove_missing)
    method <- match.arg(method, allowed_methods)

    ## computations
    if (is.null(y)) {
        out <- wdm_mat_cpp(x, method, weights, remove_missing)
        colnames(out) <- rownames(out) <- colnames(x)
    } else if (NCOL(x) == 1) {
        out <- wdm_cpp(x, y, method, weights, remove_missing)
    } else {
        out <- matrix(NA, ncol(x), ncol(y))
        for (i in seq_len(ncol(x))) {
            for (j in seq_len(ncol(y))) {
                out[i, j] <- wdm(x[, i], y[, j], method, weights, remove_missing)
            }
        }
        rownames(out) <- colnames(x)
        colnames(out) <- colnames(y)
    }

    out[is.nan(out)] <- NA
    out
}

allowed_methods <- c("pearson", "kendall", "spearman", "hoeffding", "blomqvist")

check_wdm_inputs <- function(x, y, weights, remove_missing) {
    if (!is.matrix(x) && is.null(y))
        stop("supply both 'x' and 'y' or a matrix-like 'x'")
    if (!(is.numeric(x) || is.logical(x)))
        stop("'x' must be numeric")
    if (!is.numeric(weights) | (NCOL(weights) != 1))
        stop("'weights' must be a numeric vector")
    stopifnot(is.atomic(x))
    stopifnot(is.atomic(weights))
    if (!is.null(y)) {
        if (!(is.numeric(y) || is.logical(y)))
            stop("'y' must be numeric")
        stopifnot(is.atomic(y))
        if ((NROW(x) != NROW(y)))
            stop("'x' and 'y' must have the same number of rows")
    }
    if (!is.logical(remove_missing))
        stop("remove_missing must be logical.")
}
