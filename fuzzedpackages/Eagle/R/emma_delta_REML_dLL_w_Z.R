emma.delta.REML.dLL.w.Z <- function (logdelta, lambda, etas.1, n, t1, etas.2.sq)
{
    t <- t1
    tq <- length(etas.1)
    nq <- n - t + tq
    delta <- exp(logdelta)
    etasq <- etas.1 * etas.1
    ldelta <- lambda + delta
    return(0.5 * (nq * (sum(etasq/(ldelta * ldelta)) + etas.2.sq/(delta *
        delta))/(sum(etasq/ldelta) + etas.2.sq/delta) - (sum(1/ldelta) +
        (n - t)/delta)))
}


