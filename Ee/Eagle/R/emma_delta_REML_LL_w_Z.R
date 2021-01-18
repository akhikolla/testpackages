emma.delta.REML.LL.w.Z <- function (logdelta, lambda, etas.1, n, t, etas.2.sq)
{
    tq <- length(etas.1)
    nq <- n - t + tq
    delta <- exp(logdelta)
    return(0.5 * (nq * (log(nq/(2 * pi)) - 1 - log(sum(etas.1 *
        etas.1/(lambda + delta)) + etas.2.sq/delta)) - (sum(log(lambda +
        delta)) + (n - t) * logdelta)))
}


