#' Distance between rows of two matrices
#'
#' Calculate the distances between the rows of one matrix and the rows
#' of a second matrix.
#'
#' @param x A numeric matrix.
#' @param y A second numeric matrix, with the same number of columns as `x`.
#' @param distance Indicates whether to use Euclidean distance
#'     (`"rmsd"` for root mean square difference) or the mean absolute
#'     difference (`"mad"`).
#' @param align_cols If TRUE, align the columns in the two matrices by
#' the column names.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If `0`, use [parallel::detectCores()].)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by [parallel::makeCluster()].
#'
#' @return If `x` is `m` by `p` and `y` is `n` by `p`, then the result
#' is an `m` by `n` matrix whose (i,j)th element is the distance
#' between the ith row of `x` and the jth row of `y`.
#'
#' @examples
#' p <- 10
#' n <- 5
#' m <- 3
#' x <- matrix(stats::rnorm(n*p), ncol=p)
#' rownames(x) <- LETTERS[1:n]
#' y <- matrix(stats::rnorm(m*p), ncol=p)
#' rownames(y) <- letters[1:m]
#'
#' d <- dist_betw_matrices(x, y)
#'
#' @seealso [corr_betw_matrices()], [dist_betw_arrays()]
#' @export
dist_betw_matrices <-
    function(x,y, distance=c("rmsd", "mad"), align_cols=TRUE, cores=1)
{
    distance <- match.arg(distance)

    if(!is.matrix(x)) x <- as.matrix(x)
    if(!is.matrix(y)) y <- as.matrix(y)

    # align columns by their names
    if(align_cols) {
        aligned <- align_matrix_cols(x, y)
        x <- aligned$x
        y <- aligned$y
        if(ncol(x) < 2) {
            stop("In trying to align columns, we omitted all but 1 column")
        }
    }

    if(ncol(x) != ncol(y))
        stop("x and y should have the same number of columns")

    # the C++ functions need the transpose of x and y
    x <- t(x)
    y <- t(y)

    # set up cluster
    cores <- setup_cluster(cores)

    if(distance=="rmsd") {
        func <- rmsd_betw_matrices
    } else {
        func <- mad_betw_matrices
    }

    if(n_cores(cores)==1) {
        result <- func(x, y)
    }
    else { # multi-core
        func_by_xcol <- function(i) func(x[,i,drop=FALSE], y)

        result <- cluster_lapply(cores, 1:ncol(x), func_by_xcol)
        result <- matrix(unlist(result), ncol=ncol(y), byrow=TRUE)
    }

    dimnames(result) <- list(colnames(x), colnames(y))

    result
}
