
emma.MLE <- function (y, X, K, Z = NULL, ngrids = 100, llim = -10, ulim = 10,
    esp = 1e-10, eig.L = NULL, eig.R = NULL )
{
    n <- length(y)
    t <- nrow(K)
    q <- ncol(X)
    stopifnot(ncol(K) == t)
    stopifnot(nrow(X) == n)
    if (det(crossprod(X, X)) == 0) {
        warning("X is singular")
        return(list(ML = 0, delta = 0, ve = 0, vg = 0))
    }
    if (is.null(Z)) {
        if (is.null(eig.L)) {
            eig.L <- emma.eigen.L.wo.Z(K )
        }
        if (is.null(eig.R)) {
            eig.R <- emma.eigen.R.wo.Z(K, X )
        }
        etas <- crossprod(eig.R$vectors, y)
        logdelta <- (0:ngrids)/ngrids * (ulim - llim) + llim
        m <- length(logdelta)
        delta <- exp(logdelta)
        Lambdas <- matrix(eig.R$values, n - q, m) + matrix(delta,
            n - q, m, byrow = TRUE)
        Xis <- matrix(eig.L$values, n, m) + matrix(delta, n,
            m, byrow = TRUE)
        Etasq <- matrix(etas * etas, n - q, m)
        LL <- 0.5 * (n * (log(n/(2 * pi)) - 1 - log(colSums(Etasq/Lambdas))) -
            colSums(log(Xis)))
        dLL <- 0.5 * delta * (n * colSums(Etasq/(Lambdas * Lambdas))/colSums(Etasq/Lambdas) -
            colSums(1/Xis))
        optlogdelta <- vector(length = 0)
        optLL <- vector(length = 0)
        if (dLL[1] < esp) {
            optlogdelta <- append(optlogdelta, llim)
            optLL <- append(optLL, emma.delta.ML.LL.wo.Z(llim,
                eig.R$values, etas, eig.L$values))
        }
        if (dLL[m - 1] > 0 - esp) {
            optlogdelta <- append(optlogdelta, ulim)
            optLL <- append(optLL, emma.delta.ML.LL.wo.Z(ulim,
                eig.R$values, etas, eig.L$values))
        }
        for (i in 1:(m - 1)) {
            if ((dLL[i] * dLL[i + 1] < 0 - esp * esp) && (dLL[i] >
                0) && (dLL[i + 1] < 0)) {
                r <- uniroot(emma.delta.ML.dLL.wo.Z, lower = logdelta[i],
                  upper = logdelta[i + 1], lambda = eig.R$values,
                  etas = etas, xi = eig.L$values)
                optlogdelta <- append(optlogdelta, r$root)
                optLL <- append(optLL, emma.delta.ML.LL.wo.Z(r$root,
                  eig.R$values, etas, eig.L$values))
            }
        }
    }
    else {

        if (is.null(eig.L)) {
            eig.L <- emma.eigen.L.w.Z(Z, K )
        }
        if (is.null(eig.R)) {
            eig.R <- emma.eigen.R.w.Z(Z, K, X)
        }
        etas <- crossprod(eig.R$vectors, y)
        etas.1 <- etas[1:(t - q)]
        etas.2 <- etas[(t - q + 1):(n - q)]
        etas.2.sq <- sum(etas.2 * etas.2)
        logdelta <- (0:ngrids)/ngrids * (ulim - llim) + llim
        m <- length(logdelta)
        delta <- exp(logdelta)
        Lambdas <- matrix(eig.R$values, t - q, m) + matrix(delta,
            t - q, m, byrow = TRUE)
        Xis <- matrix(eig.L$values, t, m) + matrix(delta, t,
            m, byrow = TRUE)
        Etasq <- matrix(etas.1 * etas.1, t - q, m)
        dLL <- 0.5 * delta * (n * (colSums(Etasq/(Lambdas * Lambdas)) +
            etas.2.sq/(delta * delta))/(colSums(Etasq/Lambdas) +
            etas.2.sq/delta) - (colSums(1/Xis) + (n - t)/delta))
        optlogdelta <- vector(length = 0)
        optLL <- vector(length = 0)
        if (dLL[1] < esp) {
            optlogdelta <- append(optlogdelta, llim)
            optLL <- append(optLL, emma.delta.ML.LL.w.Z(llim,
                eig.R$values, etas.1, eig.L$values, n, etas.2.sq))
        }
        if (dLL[m - 1] > 0 - esp) {
            optlogdelta <- append(optlogdelta, ulim)
            optLL <- append(optLL, emma.delta.ML.LL.w.Z(ulim,
                eig.R$values, etas.1, eig.L$values, n, etas.2.sq))
        }
        for (i in 1:(m - 1)) {
            if ((dLL[i] * dLL[i + 1] < 0 - esp * esp) && (dLL[i] >
                0) && (dLL[i + 1] < 0)) {
                r <- uniroot(emma.delta.ML.dLL.w.Z, lower = logdelta[i],
                  upper = logdelta[i + 1], lambda = eig.R$values,
                  etas.1 = etas.1, xi.1 = eig.L$values, n = n,
                  etas.2.sq = etas.2.sq)
                optlogdelta <- append(optlogdelta, r$root)
                optLL <- append(optLL, emma.delta.ML.LL.w.Z(r$root,
                  eig.R$values, etas.1, eig.L$values, n, etas.2.sq))
            }
        }
    }
    maxdelta <- exp(optlogdelta[which.max(optLL)])
    maxLL <- max(optLL)
    if (is.null(Z)) {
        maxva <- sum(etas * etas/(eig.R$values + maxdelta))/n
    }
    else {
        maxva <- (sum(etas.1 * etas.1/(eig.R$values + maxdelta)) +
            etas.2.sq/maxdelta)/n
    }
    maxve <- maxva * maxdelta
    return(list(ML = maxLL, delta = maxdelta, ve = maxve, vg = maxva))
}



