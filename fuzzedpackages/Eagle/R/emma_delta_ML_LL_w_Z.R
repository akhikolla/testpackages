
emma.delta.ML.LL.w.Z <-  function (logdelta, lambda, etas.1, xi.1, n, etas.2.sq)
{
    t <- length(xi.1)
    delta <- exp(logdelta)
    return(0.5 * (n * (log(n/(2 * pi)) - 1 - log(sum(etas.1 *
        etas.1/(lambda + delta)) + etas.2.sq/delta)) - (sum(log(xi.1 +
        delta)) + (n - t) * logdelta)))
}


