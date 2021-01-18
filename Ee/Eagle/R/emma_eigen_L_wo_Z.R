emma.eigen.L.wo.Z <- function (K)
{

    eig <- eigen(K, symmetric = TRUE)


    return(list(values = eig$values, vectors = eig$vectors))
}


