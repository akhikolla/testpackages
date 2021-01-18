#' Compute the PCA projection
#'
#' @param C data matrix used for PCA projection
#' @param L number for the top principal components
#' @import irlba irlba
#' @importFrom stats qnorm
#' @export
pca_projection_R <- function(C, L) {
    if (L >= min(dim(C))){
        eigen_res <- eigen(C)

        U <- eigen_res$vector
        V <- eigen_res$value
        eig_sort <- sort(V, decreasing = T, index.return = T)
        eig_idx <- eig_sort$ix

        W <- U[, eig_idx[1:L]]
        return (W)
    } else{
        initial_v <- as.matrix(qnorm(1:(ncol(C) + 1)/(ncol(C) + 1))[1:ncol(C)])
        eigen_res <- irlba::irlba(C, nv = L, v = initial_v)
        U <- eigen_res$u
        V <- eigen_res$v
        return (V)
    }
}

#' Get the top L eigenvalues
#' @param C data matrix used for eigendecomposition
#' @param L number for the top eigenvalues
#' @import irlba irlba
#' @export
get_major_eigenvalue <- function(C, L) {
    if (L >= min(dim(C))){
        return (base::norm(C, '2')^2);
    }else{
        #message("using irlba")
        initial_v <- as.matrix(qnorm(1:(ncol(C) + 1)/(ncol(C) + 1))[1:ncol(C)])
        eigen_res <- irlba(C, nv = L, v = initial_v)
        return (max(abs(eigen_res$v)))
    }
    #     eig_sort <- sort(V, decreasing = T, index.return = T)
    #     eig_idx <- eig_sort$ix
    #
    #     W <- U[, eig_idx[1:L]]
}

#' calculate the square distance between a, b
#' @param a a matrix with \eqn{D \times N} dimension
#' @param b a matrix with \eqn{D \times N} dimension
#' @return a numeric value for the different between a and b
#' @export
sqdist_R <- function(a, b) {
    aa <- colSums(a^2)
    bb <- colSums(b^2)
    ab <- t(a) %*% b

    aa_repmat <- matrix(rep(aa, times = ncol(b)), ncol = ncol(b), byrow = F)
    bb_repmat <- matrix(rep(bb, times = ncol(a)), nrow = ncol(a), byrow = T)
    dist <- abs(aa_repmat + bb_repmat - 2 * ab)
}

#' Perform DDRTree construction
#' @param X a matrix with \eqn{\mathbf{D \times N}} dimension which is needed to perform DDRTree construction
#' @param initial_method a function to take the data transpose of X as input and then output the reduced dimension,
#' row number should not larger than observation and column number should not be larger than variables (like isomap may only
#' return matrix on valid sample sets). Sample names of returned reduced dimension should be preserved.
#' @param dimensions reduced dimension
#' @param maxIter maximum iterations
#' @param sigma bandwidth parameter
#' @param lambda regularization parameter for inverse graph embedding
#' @param ncenter number of nodes allowed in the regularization graph
#' @param param.gamma regularization parameter for k-means (the prefix of 'param' is used to avoid name collision with gamma)
#' @param tol relative objective difference
#' @param verbose emit extensive debug output
#' @param ... additional arguments passed to DDRTree
#' @importFrom stats kmeans
#' @return a list with W, Z, stree, Y, history
#' W is the orthogonal set of d (dimensions) linear basis vector
#' Z is the reduced dimension space
#' stree is the smooth tree graph embedded in the low dimension space
#' Y represents latent points as the center of Z
#' @examples
#' data('iris')
#' subset_iris_mat <- as.matrix(t(iris[c(1, 2, 52, 103), 1:4])) #subset the data
#' #run DDRTree with ncenters equal to species number
#' DDRTree_res <- DDRTree(subset_iris_mat, dimensions = 2, maxIter = 5, sigma = 1e-2,
#' lambda = 1, ncenter = 3, param.gamma = 10, tol = 1e-2, verbose = FALSE)
#' Z <- DDRTree_res$Z #obatain matrix
#' Y <- DDRTree_res$Y
#' stree <- DDRTree_res$stree
#' plot(Z[1, ], Z[2, ], col = iris[c(1, 2, 52, 103), 'Species']) #reduced dimension
#' legend("center", legend = unique(iris[c(1, 2, 52, 103), 'Species']), cex=0.8,
#' col=unique(iris[c(1, 2, 52, 103), 'Species']), pch = 1) #legend
#' title(main="DDRTree reduced dimension", col.main="red", font.main=4)
#' dev.off()
#' plot(Y[1, ], Y[2, ], col = 'blue', pch = 17) #center of the Z
#' title(main="DDRTree smooth principal curves", col.main="red", font.main=4)
#'
#' #run DDRTree with ncenters equal to species number
#' DDRTree_res <- DDRTree(subset_iris_mat, dimensions = 2, maxIter = 5, sigma = 1e-3,
#' lambda = 1, ncenter = NULL,param.gamma = 10, tol = 1e-2, verbose = FALSE)
#' Z <- DDRTree_res$Z #obatain matrix
#' Y <- DDRTree_res$Y
#' stree <- DDRTree_res$stree
#' plot(Z[1, ], Z[2, ], col = iris[c(1, 2, 52, 103), 'Species']) #reduced dimension
#' legend("center", legend = unique(iris[c(1, 2, 52, 103), 'Species']), cex=0.8,
#' col=unique(iris[c(1, 2, 52, 103), 'Species']), pch = 1) #legend
#' title(main="DDRTree reduced dimension", col.main="red", font.main=4)
#' dev.off()
#' plot(Y[1, ], Y[2, ], col = 'blue', pch = 2) #center of the Z
#' title(main="DDRTree smooth principal graphs", col.main="red", font.main=4)
#' @export
#'
DDRTree <- function(X,
                        dimensions = 2,
                        initial_method = NULL,
                        maxIter = 20,
                        sigma = 1e-3,
                        lambda = NULL,
                        ncenter = NULL,
                        param.gamma = 10,
                        tol = 1e-3,
                        verbose = F, ...) {

    D <- nrow(X)
    N <- ncol(X)

    #initialization
    W <- pca_projection_R(X %*% t(X), dimensions)
    if(is.null(initial_method)){
        Z <- t(W) %*% X
    }
    else{
      tmp <- initial_method(X, ...) #a function to return reduced dimension data
      if(ncol(tmp) > D | nrow(tmp) > N)
        stop('The dimension reduction method passed need to return correct dimensions')
      Z <- tmp[, 1:dimensions]
      Z <- t(Z)
    }

    if(is.null(ncenter)) {
        K <- N
        Y <- Z[, 1:K]
    }
    else {
        K <- ncenter
        if (K > ncol(Z))
            stop("Error: ncenters must be greater than or equal to ncol(X)")
        centers = t(Z)[seq(1, ncol(Z), length.out=K),]
        kmean_res <- kmeans(t(Z), K, centers=centers)
        Y <- kmean_res$centers
        Y <- t(Y)
    }

    if (is.null(lambda)){
        lambda = 5 * ncol(X)
    }
    ddrtree_res <- DDRTree_reduce_dim(X, Z, Y, W, dimensions, maxIter, K,  sigma,  lambda,  param.gamma, tol, verbose)

    return(list(W = ddrtree_res$W, Z = ddrtree_res$Z, stree = ddrtree_res$stree, Y = ddrtree_res$Y, X = ddrtree_res$X, objective_vals = ddrtree_res$objective_vals, history = NULL))
}
