emma.eigen.L.w.Z <- function (Z, K, complete = TRUE )
{
    if (complete == FALSE) {
        vids <- colSums(Z) > 0
        Z <- Z[, vids]
        K <- K[vids, vids]
    }
    res <- K %*% crossprod(Z, Z)
    ## cannot use eigen_mgpu here because matrix is not symmetric
    eig <- eigen(res, symmetric = FALSE, EISPACK = TRUE)
    return(list(values = eig$values, vectors = qr.Q(qr(Z %*%
        eig$vectors), complete = TRUE)))
}

