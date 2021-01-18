#' Expected value and variance of renewal count process
#'
#' Compute numerically expected values and variances of renewal count processes.
#'
#' \code{evCount_conv_bi} computes the expected value and variance of renewal
#' count processes for the builtin distirbutions of inter-arrival times.
#' 
#' @param xmax unsigned integer maximum count to be used.
#' @param dist TODO
#' @param method TODO
#' @param distPars TODO
#' @inheritParams dCount_conv_bi
#' @return
#'     a named list with components \code{"ExpectedValue"} and \code{"Variance"}.
#' @examples
#' pwei_user <- function(tt, distP) {
#'     alpha <- exp(-log(distP[["scale"]]) / distP[["shape"]])
#'     pweibull(q = tt, scale = alpha, shape = distP[["shape"]],
#'              lower.tail = FALSE)
#' }
#'
#' ## ev convolution Poisson count
#' lambda <- 2.56
#' beta <- 1
#' distPars <- list(scale = lambda, shape = beta)
#'
#' evbi <- evCount_conv_bi(20, distPars, dist = "weibull")
#' evu <- evCount_conv_user(20, distPars, c(2, 2), pwei_user, "dePril")
#'
#' c(evbi[["ExpectedValue"]], lambda)
#' c(evu[["ExpectedValue"]], lambda )
#' c(evbi[["Variance"]], lambda     )
#' c(evu[["Variance"]], lambda      )
#'
#' ## ev convolution weibull count
#' lambda <- 2.56
#' beta <- 1.35
#' distPars <- list(scale = lambda, shape = beta)
#'
#' evbi <- evCount_conv_bi(20, distPars, dist = "weibull")
#' evu <- evCount_conv_user(20, distPars, c(2.35, 2), pwei_user, "dePril")
#'
#' x <- 1:20
#' px <- dCount_conv_bi(x, distPars, "weibull", "dePril",
#'                      nsteps = 100)
#' ev <- sum(x * px)
#' var <- sum(x^2 * px) - ev^2
#'
#' c(evbi[["ExpectedValue"]], ev)
#' c(evu[["ExpectedValue"]], ev )
#' c(evbi[["Variance"]], var    )
#' c(evu[["Variance"]], var     )
#' @export
evCount_conv_bi <- function(xmax, distPars,
                            dist = c("weibull", "gamma", "gengamma",
                                "burr"),
                            method = c( "dePril", "direct", "naive"),
                            nsteps = 100, time = 1.0,
                            extrap = TRUE) {

    dist <- match.arg(dist)
    method <- match.arg(method)

    x <- 1:xmax
    px <- dCount_conv_bi(x, distPars, dist, method, nsteps, time,
                         extrap, log = FALSE)

    ## adjust probability 
    px[px < 0] <- 0
        
    ev <- sum(x * px)
    ev2 <- sum(x^2 * px)
    var <- ev2 - ev^2
    if (var < 0) {
        px <- zapsmall(px)
        ##px <- px / sum(px)
        ev <- sum(x * px)
        ev2 <- sum(x^2 * px)
        var <- ev2 - ev^2
    }
    
    list(ExpectedValue = ev, Variance = var)
}

#' % renewal count expected value and variance (user)
#'
#' % Expected value and variance of the renewal count process computed numerically
#'
#' \code{evCount_conv_user} computes the expected value and variance for a user
#' specified distirbution of the inter-arrival times.
#' 
#' @param extrapolPars ma::vec of length 2. The extrapolation values.
#' @param survR function, user supplied survival function; should have signature
#'     \code{function(t, distPars)}, where \code{t} is a positive real number
#'     (the time where the survival function is evaluated) and \code{distPars}
#'     is a list of distribution parameters. It should return a double value.
#' @inheritParams evCount_conv_bi
#' @return % List named vector with ExpectedValue and Variance
#' @rdname evCount_conv_bi
#' @export
evCount_conv_user <- function(xmax, distPars, extrapolPars, survR,
                              method = c( "dePril", "direct", "naive"),
                              nsteps = 100, time = 1.0,
                              extrap = TRUE) {

    method <- match.arg(method)

    x <- 1:xmax
    px <- dCount_conv_user(x, distPars, extrapolPars, survR,
                           method, nsteps, time,
                           extrap, log = FALSE)

    ## adjust probability 
    px[px < 0] <- 0
    
    ev <- sum(x * px)
    ev2 <- sum(x^2 * px)
    var <- ev2 - ev^2
    if (var < 0) {
        px <- zapsmall(px)
        ##px <- px / sum(px)
        ev <- sum(x * px)
        ev2 <- sum(x^2 * px)
        var <- ev2 - ev^2
    }
    
    list(ExpectedValue = ev, Variance = var)
}
