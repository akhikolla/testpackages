##' Calculate entropy of Truncated Normal distributions
##'
##' @title Truncated Normal Distribution Entropy
##' @param .mean Length K vector with the means of the K Normal
##' distributions prior to truncation
##' @param .sd Length K vector with the standard deviations of the K Normal
##' distributions prior to truncation
##' @param .low Length K vector with the lower truncation bound of the K
##' Normal distributions prior to truncation
##' @param .high Length K vector with the upper truncation bound of the K
##' Normal distributions prior to truncation
##' @return Length K vector with the entropies associated with each of the K
##' Truncated Normal distributions
##' @author Jonathan Olmsted
##' @examples
##' lows <- c(-1, 5, -100, 4, 4, -100, 7)
##' highs <- c(1, 100, 10, 7, 4.1, 100, 100)
##' enttn(.mean = rep(0, length(lows)),
##'       .sd = rep(1, length(lows)),
##'       .low = lows,
##'       .high = highs
##'       )
##'

enttn <- function(.mean = rep(0, 1),
                  .sd = rep(1, length(.mean)),
                  .low = rep(-Inf, length(.mean)),
                  .high = rep(Inf, length(.mean))
                  ) {
    out <- .Call("enttnRcpp",
                 mean = .mean,
                 sd = .sd,
                 low = .low,
                 high = .high
                 )
    return(out)
}
