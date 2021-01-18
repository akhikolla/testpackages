#' Uniform scaled beta distribution
#'
#' The uniform scaled beta (USB) distribution describes the distribution of
#' the random variable
#' \deqn{U_{b, \nu} = U + \theta(B - 0.5),}
#' where \eqn{U} is a \eqn{U[-0.5, 0.5]} random variable, \eqn{B} is a
#' \eqn{Beta(\nu, \nu)} random variable, and \eqn{theta > 0, \nu >= 1}.
#'
#' @param x vector of quantiles.
#' @param theta scale parameter of the USB distribution.
#' @param nu smoothness parameter of the USB distribution.
#'
#' @references
#' Nagler, T. (2017). *A generic approach to nonparametric function
#' estimation with mixed data.* [arXiv:1704.07457](https://arxiv.org/pdf/1704.07457.pdf)
#'
#' @examples
#' # plot distribution
#' sq <- seq(-0.8, 0.8, by = 0.01)
#' plot(sq, dusb(sq), type = "l")
#' lines(sq, dusb(sq, theta = 0.25), col = 2)
#' lines(sq, dusb(sq, theta = 0.25, nu = 10), col = 3)
#'
#' # simulate from the distribution
#' x <- rusb(100, theta = 0.3, nu = 0)
#'
#' @importFrom stats pbeta
#' @export
dusb <- function(x, theta = 0, nu = 5) {
    stopifnot(theta >= 0)
    stopifnot(nu > 0)
    if (theta == 0)
        return(as.numeric(abs(x) < 0.5))
    pbeta((x + 0.5) / theta + 0.5, nu, nu) - pbeta((x - 0.5) / theta + 0.5, nu, nu)
}

#' @rdname dusb
#' @param n number of observations.
#' @param quasi logical indicating whether quasi random numbers
#'   ([qrng::ghalton()]) should be used for generating uniforms (which are then
#'   transformed by the quantile function)
#' @importFrom qrng ghalton
#' @importFrom stats qbeta rbeta
#' @export
rusb <- function(n, theta = 0, nu = 5, quasi = FALSE) {
    stopifnot(theta >= 0)
    stopifnot(theta <= 1)
    if (!quasi) {
        x <- (runif(n) - 0.5)
        if (theta > 0)
            x <- x + theta * (rbeta(n, nu, nu) - 0.5)
    } else {
        # permute the quasi-random sequence randomly to avoid correlation
        # between several usb random variables
        if (theta == 0) {
            x <- (qrng::ghalton(n, d = 1) - 0.5)[sample(n)]
        } else {
            u <- qrng::ghalton(n, d = 2)[sample(n), ]
            x <- (u[, 1]  - 0.5) + theta * (qbeta(u[, 2], nu, nu) - 0.5)
        }
    }

    x
}
