#' Continuous convolution density estimator
#'
#' The continuous convolution kernel density estimator is defined as the
#' classical kernel density estimator based on continuously convoluted data (see
#' [cont_conv()]). [cckde()] fits the estimator (including bandwidth selection),
#' [dcckde()] and [predict.cckde()] can be used to evaluate the estimator.
#'
#' @param x a matrix or data frame containing the data (or evaluation points).
#' @param bw vector of bandwidth parameter; if `NULL`, the bandwidths are
#'   selected automatically by likelihood cross validation.
#' @param mult bandwidth multiplier; either a positive number or a vector of
#'   such. Each bandwidth parameter is multiplied with the corresponding
#'   multiplier.
#' @param theta scale parameter of the USB distribution (see, [dusb()]).
#' @param nu smoothness parameter of the USB distribution (see, [dusb()]).
#'   The estimator uses the Epanechnikov kernel for smoothing and the USB
#'   distribution for continuous convolution (default parameters correspond to
#'   the uniform distribution on \eqn{[-0.5, 0.5]}.
#' @param object `cckde` object.
#' @param newdata matrix or data frame containing evaluation points.
#' @param ... unused.
#'
#' @details If a variable should be treated as ordered discrete, declare it as
#'   [ordered()], factors are expanded into discrete dummy codings.
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
#' fit <- cckde(dat)  # fit estimator
#' dcckde(dat, fit)   # evaluate density
#' predict(fit, dat)  # equivalent
#'
#' @export
#' @useDynLib cctools
cckde <- function(x, bw = NULL, mult = 1, theta = 0, nu = 5, ...) {
    # continuous convolution of the data
    x_cc <- cont_conv(x, theta = theta, nu = nu)

    if (is.null(bw)) {
        # find optimal bandwidths using likelihood cross-validation
        x_eval <- expand_as_numeric(x)
        bw <- select_bw(x = x_eval,
                        x_cc = x_cc,
                        i_disc = attr(x_cc, "i_disc"),
                        bw_min = 0.5 - theta)
    }

    # adjust bws
    stopifnot(all(bw > 0))
    stopifnot(all(mult > 0))
    bw <- expand_vec(bw, x) * expand_vec(mult, x)
    names(bw) <- colnames(x_cc)

    # create and return cckde object
    structure(
        list(x_cc = x_cc, bw = bw, theta = theta, ell = nu),
        class = "cckde"
    )
}


#' @rdname cckde
#' @export
dcckde <- function(x, object) {
    stopifnot(inherits(object, "cckde"))
    # must be numeric, factors are expanded
    x <- expand_as_numeric(x)
    # variables must be in same order
    x <- x[, colnames(object$x_cc), drop = FALSE]

    # raw density
    f <- c(eval_mvkde(x, object$x_cc, object$bw))

    ## normalize such that discrete variables' marginal densities sum to 1
    i_disc <- attr(object$x_cc, "i_disc")
    if (length(i_disc) > 0) {
        for (i in i_disc) {
            lvls <- seq_along(attr(object$x_cc, "levels")[[i]]) - 1
            f <- f /
                sum(eval_mvkde(as.matrix(lvls),
                               object$x_cc[, i, drop = FALSE],
                               object$bw[i, drop = FALSE]))
        }
    }

    f
}

#' @rdname cckde
#' @export
predict.cckde <- function(object, newdata, ...)
    dcckde(newdata, object)

#' @importFrom stats IQR optim pbeta rbeta runif sd
#' @importFrom Rcpp evalCpp
#' @noRd
select_bw <- function(x, x_cc, i_disc = integer(0), bw_min = 0) {
    ## set lower bounds for the bandwidth of each variable
    bw_lower <- rep(0, ncol(x))
    bw_lower[i_disc] <- bw_min

    ## set starting values by normal reference rule
    # - 2.3449 is the constant for Epanechnikov kernel
    # - n^(-1 / (4 + d_cont)) is the optimal order, where d_cont is the number
    #   of continuous variables
    # - we use d_cont / 2 for the middle way between asymptotics and small
    #   sample behavior
    bw_start_fun <- function(y)
        2.3449 * sd(y) * nrow(x_cc)^(-1 / (4 + ncol(x_cc) - length(i_disc) / 2))
    bw_start <- apply(x_cc, 2, bw_start_fun)
    bw_start <- pmax(bw_start, bw_lower)  # adjust with lower bounds

    ## find optimal bandwidth by likelihood cross-validation
    opt <- optim(
        bw_start,
        function(bw) lcv_mvkde_disc(x, x_cc, bw),
        lower = bw_lower,
        method = "L-BFGS-B",
        control = list(fnscale = -1, parscale = apply(x_cc, 2, sd), pgtol = 0)
    )

    ## return optimal bandwidths
    opt$par
}
