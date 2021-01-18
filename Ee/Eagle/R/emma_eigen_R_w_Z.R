
emma.eigen.R.w.Z <-  function (Z, K, X, complete = TRUE  )
{


    if (complete == FALSE) {
        vids <- colSums(Z) > 0
        Z <- Z[, vids]
        K <- K[vids, vids]
    }
    n <- nrow(Z)
    t <- ncol(Z)
    q <- ncol(X)
    # matrix too large for GPU based code
     #SZ <- Z - X %*%   solve(crossprod(X, X))    %*% crossprod(X, Z)
     SZ <- Z - X %*%    chol2inv(chol(crossprod(X, X)))    %*% crossprod(X, Z)

    eig <- eigen(K %*% crossprod(Z, SZ), symmetric = FALSE, EISPACK = TRUE)
    



    if (is.complex(eig$values)) {
        eig$values <- Re(eig$values)
        eig$vectors <- Re(eig$vectors)
    }

    qr.X <- qr.Q(qr(X))
    values = eig$values[1:(t - q)]
    vectors = qr.Q(qr(cbind(SZ %*% eig$vectors[, 1:(t - q)], qr.X)), complete = TRUE)[, c(1:(t - q), (t + 1):n)]
   return(list(values=values, vectors=vectors))

}


