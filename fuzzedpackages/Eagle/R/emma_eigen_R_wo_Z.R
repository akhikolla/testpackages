
emma.eigen.R.wo.Z <-  function (K, X )
{
    n <- nrow(X)
    q <- ncol(X)
    # S <- -1 *  X %*% solve(crossprod(X, X)) %*% t(X)
     S <- -1 *  X %*%  chol2inv(chol(crossprod(X, X)))  %*% t(X)

    diag(S) <- diag(S) + 1
    diag(K) <- diag(K) + 1  ## old code: K + dn
    gc()

       eig <- eigen(S %*% K %*% S, symmetric = TRUE)


    stopifnot(!is.complex(eig$values))
    return(list(values = eig$values[1:(n - q)] - 1, vectors = eig$vectors[,
        1:(n - q)]))
}


