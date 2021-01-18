##' Calculate density of Truncated Normal distributions
##'
##' @title Truncated Normal Distribution Density
##' @param .x Length K vector of the points at which to evaluate the density
##' @param .mean Length K vector with the means of the K Normal distributions
##'     *prior* to truncation
##' @param .sd Length K vector with the standard deviations of the K Normal
##'     distributions *prior* to truncation
##' @param .low Length K vector with the lower truncation bound of the K Normal
##'     distributions *prior* to truncation
##' @param .high Length K vector with the upper truncation bound of the K Normal
##'     distributions *prior* to truncation
##' @param .checks Logical indicating whether inputs and outputs should be
##'     checked and either stop (for bad inputs) or warn (for likely bad outputs)
##' @return Length K vector with the entropies associated with each of the K
##'     Truncated Normal distributions
##' @author Jonathan Olmsted
##' @examples
##' lows <- c(-1, 5, -100, 4, 4, -100, 7)
##' highs <- c(1, 100, 10, 7, 4.1, 100, 100)
##' dtn(.x = rep(0, length(lows)),
##'     .mean = rep(0, length(lows)),
##'     .sd = rep(1, length(lows)),
##      .low = lows,
##'     .high = highs
##'     )

dtn <- function(.x = 0,
                .mean = rep(0, length(.x)),
                .sd = rep(1, length(.x)),
                .low = rep(-Inf, length(.x)),
                .high = rep(Inf, length(.x)),
                .checks = TRUE
                ) {
    if (.checks) {
        checkInputs(.mean, .sd, .low, .high)
    }
    out <- .Call("dtnRcpp",
                 x = .x,
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
