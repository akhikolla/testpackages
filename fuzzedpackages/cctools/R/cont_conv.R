#' Continuous convolution
#'
#' Applies the continuous convolution trick, i.e. adding continuous noise to all
#' discrete variables. If a variable should be treated as discrete, declare it
#' as [ordered()] (passed to [expand_as_numeric()]).
#'
#' @param x data; numeric matrix or data frame.
#' @param theta scale parameter of the USB distribution (see, [dusb()]).
#' @param nu smoothness parameter of the USB distribution (see, [dusb()]). The
#'   estimator uses the Epanechnikov kernel for smoothing and the USB for
#'   continuous convolution (default parameters correspond to the \eqn{U[-0.5,
#'   0.5]} distribution).
#' @param quasi logical indicating whether quasi random numbers sholuld be used
#'   ([qrng::ghalton()]); only works for `theta = 0`.
#'
#' @return A data frame with noise added to each discrete variable (ordered
#'   columns).
#'
#' @details The UPSB distribution ([dusb()]) is used as the noise distribution.
#'   Discrete variables are assumed to be integer-valued.
#'
#' @references
#' Nagler, T. (2017). *A generic approach to nonparametric function
#' estimation with mixed data.* [arXiv:1704.07457](https://arxiv.org/pdf/1704.07457.pdf)
#'
#' @examples
#' # dummy data with discrete variables
#' dat <- data.frame(
#'     F1 = factor(rbinom(10, 4, 0.1), 0:4),
#'     Z1 = ordered(rbinom(10, 5, 0.5), 0:5),
#'     Z2 = ordered(rpois(10, 1), 0:10),
#'     X1 = rnorm(10),
#'     X2 = rexp(10)
#' )
#'
#' pairs(dat)
#' pairs(expand_as_numeric(dat))  # expanded variables without noise
#' pairs(cont_conv(dat))          # continuously convoluted data
#'
#' @export
cont_conv <- function(x, theta = 0, nu = 5, quasi = TRUE) {
    if (inherits(x, "expanded_as_numeric"))
        return(x)
    cc_add_noise(expand_as_numeric(x), theta, nu, quasi)
}

#' Add noise to discrete variables
#'
#' @param x data matrix created by `expand_as_numeric()` (this is important,
#'   because we need to know which variables are discrete).
#' @noRd
cc_add_noise <- function(x, theta, nu, quasi) {
    stopifnot(inherits(x, "expanded_as_numeric"))
    i_disc <- attr(x, "i_disc")
    n_disc <- length(i_disc)
    if (n_disc > 0) {
        E <- matrix(rusb(n_disc * nrow(x), theta, nu, quasi), nrow(x), n_disc)
        x[, i_disc] <- x[, i_disc] + E
    }

    attr(x, "theta") <- theta
    attr(x, "nu") <- nu
    x
}

