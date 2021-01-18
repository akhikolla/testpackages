emma.delta.ML.dLL.w.Z <-  function (logdelta, lambda, etas.1, xi.1, n, etas.2.sq)
{
    t <- length(xi.1)
    delta <- exp(logdelta)
    etasq <- etas.1 * etas.1
    ldelta <- lambda + delta
    return(0.5 * (n * (sum(etasq/(ldelta * ldelta)) + etas.2.sq/(delta *
        delta))/(sum(etasq/ldelta) + etas.2.sq/delta) - (sum(1/(xi.1 +
        delta)) + (n - t)/delta)))
}


