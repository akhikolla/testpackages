#' Distance between rows of two arrays
#'
#' Calculate the distances between the rows of two multi-dimensional
#' arrays.
#'
#' @param x A numeric array.
#' @param y A second numeric array, with the same dimensions as `x`.
#' @param distance Indicates whether to use Euclidean distance
#'     (`"rmsd"` for root mean square difference) or the mean absolute
#'     difference (`"mad"`).
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If `0`, use [parallel::detectCores()].)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by [parallel::makeCluster()].
#'
#' @details
#' The two arrays need to have the same dimensions, except for the
#' leading dimension (rows). They are turned into matrices by merging
#' all but the leading dimension, and then they're sent to
#' [dist_betw_matrices()].
#'
#' @return If `x` and `y` have `m` and `n` rows, respectively, the
#' result is an `m` by `n` matrix whose (i,j)th element is the
#' distance between the ith column of `x` and the jth column of
#' `y`.
#'
#' @examples
#' p <- 10
#' k <- 6
#' n <- 5
#' m <- 3
#' x <- array(stats::rnorm(n*k*p), dim=c(n,k,p))
#' rownames(x) <- LETTERS[1:n]
#' y <- array(stats::rnorm(m*k*p), dim=c(m,k,p))
#' rownames(y) <- letters[1:m]
#'
#' d <- dist_betw_arrays(x, y)
#'
#' @seealso [dist_betw_matrices()], [corr_betw_matrices()]
#' @export
dist_betw_arrays <-
    function(x,y, distance=c("rmsd", "mad"), cores=1)
{
    distance <- match.arg(distance)

    if(is.matrix(x) && is.matrix(y)) # matrices ... use dist_betw_matrices()
        return(dist_betw_matrices(x,y, distance=distance, align_cols=FALSE, cores=cores))

    if(!is.array(x) || !is.array(y))
        stop("x and y must both be arrays")

    dx <- dim(x)
    dy <- dim(y)
    if(length(dx) != length(dy))
        stop("length(dim(x)) != length(dim(y))")
    if(!all(dx[-1] == dy[-1]))
        stop("Need all but leading dimensions to match")

    # turn into matrices, preserving row names
    rnx <- rownames(x)
    rny <- rownames(y)
    dim(x) <- c(dx[1], prod(dx[-1]))
    dim(y) <- c(dy[1], prod(dy[-1]))
    rownames(x) <- rnx
    rownames(y) <- rny

    dist_betw_matrices(x, y, distance=distance, align_cols=FALSE, cores=cores)
}
