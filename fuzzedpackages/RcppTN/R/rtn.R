##' Sample from Truncated Normal distributions
##'
##' The special values of -Inf and Inf are valid values in the .low and .high
##' arguments, respectively. The implementation is from Robert (1995). The
##' computation is written in Rcpp-based C++ code, but respects R's RNG
##' state. The draws from this function are reproducible because it respects R's
##' RNG state. Draws using this algorithm (whether implemented in R code or C++)
##' will be the same if seeded correctly. However, you should not expect these
##' draws to match those from another algorithm.
##'
##' @title Truncated Normal Distribution RNG
##' @param .mean Length K vector with the means of the K Normal
##' distributions prior to truncation
##' @param .sd Length K vector with the standard deviations of the K Normal
##' distributions prior to truncation
##' @param .low Length K vector with the lower truncation bound of the K
##' Normal distributions prior to truncation
##' @param .high Length K vector with the upper truncation bound of the K
##' Normal distributions prior to truncation
##' @param .checks Length 1 logical vector indicating whether to perform checks
##' (safer) or not (faster) on the input parameters
##' @return A length K vector of expectations corresponding to the Truncated
##' Normal distributions. NAs are returned (with a warning) for invalid
##' parameter values.
##' @author Jonathan Olmsted
##' @references Robert, Christian P. ``Simulation of truncated normal
##' variables''. Statistics and Computing 5.2 (1995):
##' 121-125. \url{http://dx.doi.org/10.1007/BF00143942}
##' @examples
##' set.seed(1)
##' rtn(0, 1, -Inf, Inf) # single draw from a single distribution
##'
##' ## [1] -0.6264538
##'
##' set.seed(1)
##' rtn(0, 1, -Inf, Inf) # again, because it respects the RNG state
##'
##' ## [1] -0.6264538
##'
##' rtn(rep(0, 3),
##'     rep(1, 3),
##'     rep(-Inf, 3),
##'     rep(Inf, 3)
##'     ) # multiple draws from a single distribution
##'
##' ## [1]  0.1836433 -0.8356286  1.5952808
##'
##' rtn(c(0, 0),
##'     c(1, 1),
##'     c(-Inf, 5),
##'     c(1, Inf)
##'     ) # multiple draws, each from a different distribution
##' ## [1] 0.3295078 5.3917301
rtn <- function(.mean = rep(0, 1),
                .sd = rep(1, length(.mean)),
                .low = rep(-Inf, length(.mean)),
                .high = rep(Inf, length(.mean)),
                .checks = TRUE
                ) {
    if (.checks) {
        checkInputs(.mean, .sd, .low, .high)
    }
    out <- .Call("rtnRcpp",
                 mean = .mean,
                 sd = .sd,
                 low = .low,
                 high = .high
                 )
    if (.checks) {
        checkOutputs(out)
    }
    return(out)
}
