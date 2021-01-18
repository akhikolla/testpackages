#' @rdname degp3
#' @export
pegp3 <- function(q, kappa=1, sigma, xi, u=0, lower.tail=TRUE, log.p=FALSE){
    res <- pgpd(q, sigma=sigma, xi=xi, u=u, lower.tail=lower.tail, log.p=TRUE)
    if (isTRUE(lower.tail)) {
        res <- kappa * res
    } else {
        res <- .log1mexp(kappa * .log1mexp(res))
    }

    if (log.p) res else exp(res)
}
