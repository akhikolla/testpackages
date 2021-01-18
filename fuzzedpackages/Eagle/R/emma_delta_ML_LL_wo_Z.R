
emma.delta.ML.LL.wo.Z <- function (logdelta, lambda, etas, xi)
{
    n <- length(xi)
    delta <- exp(logdelta)
    return(0.5 * (n * (log(n/(2 * pi)) - 1 - log(sum((etas *
        etas)/(lambda + delta)))) - sum(log(xi + delta))))
}



