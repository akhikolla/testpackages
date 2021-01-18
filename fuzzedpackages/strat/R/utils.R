#' @useDynLib strat
#' @importFrom Rcpp evalCpp
NULL


wtd_rank <- function (x, weights = NULL, normwt = FALSE, na.rm = TRUE){
  if (!length(weights))
    return(rank(x, na.last = if (na.rm) NA else TRUE))
  tab <- Hmisc::wtd.table(x, weights, normwt = normwt, na.rm = na.rm)
  freqs <- tab$sum.of.weights
  r <- cumsum(freqs) - 0.5 * (freqs - 1)
  stats::approx(tab$x, r, xout = x, rule = 2)$y
}

clean <- function(outcome, strata, weights = NULL, group = NULL){

    # check input
    if (!length(weights))
        weights <- rep(1, length(outcome))
    if (!length(group))
        group <- rep(1, length(outcome))

    # error messages
    if (!is.numeric(outcome) || !length(outcome))
        stop("outcome has to be a numeric vector")
    if (!is.numeric(weights))
        stop("weights has to be a numeric vector")
    if (length(outcome) != length(strata))
        stop("outcome and strata have to be of equal length")
    if (length(outcome) != length(weights))
        stop("outcome and weights have to be of equal length")
    if (length(outcome) != length(group))
        stop("outcome and group have to be of equal length")

    # return a data frame
    ok <- stats::complete.cases(outcome, strata, weights)
    n <- sum(ok)
    if (n==0) stop("no complete cases!")
    outcome <- outcome[ok]
    strata <- factor(strata[ok])
    weights <- weights[ok]/sum(weights[ok]) * n
    group <- factor(group[ok])
    prank <- wtd_rank(outcome, weights, normwt = TRUE, na.rm=TRUE)/n
    out <- data.frame(prank, strata, weights, group)
    return(out)
}

.onUnload <- function (libpath) {
  library.dynam.unload("strat", libpath)
}
