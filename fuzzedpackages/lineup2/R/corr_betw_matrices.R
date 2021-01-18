#' Calculate correlations between columns of two matrices
#'
#' For matrices x and y, calculate the correlation between columns of x and
#' columns of y.
#'
#' Missing values (`NA`) are ignored, and we calculate the correlation
#' using all complete pairs, as in [stats::cor()] with
#' `use="pairwise.complete.obs"`.
#'
#' @param x A numeric matrix.
#' @param y A numeric matrix with the same number of rows as `x`.
#' @param what Indicates which correlations to calculate and return.  See
#' value, below.
#' @param corr_threshold Threshold on correlations if `what="bestpairs"`.
#' @param align_rows If TRUE, align the rows in the two matrices by
#' the row names.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If `0`, use [parallel::detectCores()].)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by [parallel::makeCluster()].
#'
#' @return If `what="paired"`, the return value is a vector of
#' correlations, between columns of `x` and the corresponding column of
#' `y`.  `x` and `y` must have the same number of columns.
#'
#' If `what="bestright"`, we return a data frame of size `ncol(x)` by
#' `3`, with the \eqn{i}th row being the maximum correlation between
#' column \eqn{i} of `x` and a column of `y`, and then the
#' `y`-column index and `y`-column name with that correlation.  (In
#' case of ties, we give the first one.)
#'
#' If `what="bestpairs"`, we return a data frame with five columns,
#' containing all pairs of columns (with one in `x` and one in `y`)
#' with correlation \eqn{\ge} `corr_threshold`.  Each row corresponds to a
#' column pair, and contains the correlation and then the `x`- and
#' `y`-column indices followed by the `x`- and `y`-column names.
#'
#' If `what="all"`, the output is a matrix of size `ncol(x)` by
#' `ncol(y)`, with all correlations between columns of `x` and
#' columns of `y`.
#'
#' @examples
#' # use the provided data, and first align the rows
#' aligned <- align_matrix_rows(lineup2ex$gastroc, lineup2ex$islet)
#'
#' # correlations for each column in x with each in y
#' result_pairs <- corr_betw_matrices(aligned[[1]], aligned[[2]], "paired")
#'
#' # subset columns to those with correlation > 0.75
#' gastroc <- lineup2ex$gastroc[,result_pairs > 0.75]
#' islet <- lineup2ex$islet[,result_pairs > 0.75]
#'
#' # similarity matrix for the two sets of rows
#' # (by transposing and using what="all")
#' corr_betw_samples <- corr_betw_matrices(t(gastroc), t(islet), "all")
#'
#' # for each column in x, find most correlated column in y
#' # (max in each row of result_all)
#' bestright <- corr_betw_matrices(t(gastroc), t(islet), "bestright")
#'
#' # correlations that exceed a threshold
#' bestpairs <- corr_betw_matrices(t(gastroc), t(islet), "bestpairs", corr_threshold=0.8)
#'
#' @seealso [dist_betw_matrices()], [dist_betw_arrays()]
#' @export
corr_betw_matrices <-
    function(x, y, what=c("paired", "bestright", "bestpairs", "all"),
             corr_threshold=0.9, align_rows=TRUE, cores=1)
{
    if(!is.matrix(x)) x <- as.matrix(x)
    if(!is.matrix(y)) y <- as.matrix(y)

    # align rows by their names
    if(align_rows) {
        aligned <- align_matrix_rows(x, y)
        x <- aligned$x
        y <- aligned$y
        if(nrow(x) < 2) {
            stop("In trying to align rows, we omitted all but 1 row")
        }
    }

    if(nrow(x) != nrow(y))
        stop("x and y should have the same number of rows")

    if(is.null(colnames(x))) colnames(x) <- paste("V", 1:ncol(x), sep="")
    if(is.null(colnames(y))) colnames(y) <- paste("V", 1:ncol(y), sep="")

    what <- match.arg(what)

    # set up cluster
    cores <- setup_cluster(cores)

    switch(what,
           paired=corr_betw_matrices_paired(x, y, cores=cores),
           bestright=corr_betw_matrices_unpaired_bestright(x,y, cores=cores),
           bestpairs=corr_betw_matrices_unpaired_bestpairs(x,y, corr_threshold, cores=cores),
           all=corr_betw_matrices_unpaired_all(x,y, cores=cores))
}

corr_betw_matrices_paired <-
    function(x,y, cores=1)
{
    px <- ncol(x)
    py <- ncol(y)

    if(py != px) stop('what="paired", but ncol(x) != ncol(y) [', px, ' != ', py, ']')

    if(n_cores(cores)==1) {
        result <- .corr_betw_matrices_paired(x, y)
    } else {
        f <- function(i) .corr_betw_matrices_paired(x[,i,drop=FALSE], y[,i,drop=FALSE])
        result <- unlist(cluster_lapply(cores, 1:ncol(x), f))
    }

    names(result) <- colnames(x)
    result
}

corr_betw_matrices_unpaired_bestright <-
    function(x,y, cores=1)
{
    if(n_cores(cores)==1) {
        result <- .corr_betw_matrices_unpaired_bestright(x, y)
        result <- as.data.frame(result)
    } else {
        f <- function(i) as.data.frame(.corr_betw_matrices_unpaired_bestright(x[,i,drop=FALSE], y))
        result <- cluster_lapply(cores, 1:ncol(x), f)
        result <- do.call("rbind", result)
    }

    colnames(result) <- c("corr", "yindex")
    result$yindex <- as.integer(result$yindex)
    rownames(result) <- colnames(x)

    cbind(result, ycol=colnames(y)[result[,2]], stringsAsFactors=FALSE)
}

corr_betw_matrices_unpaired_bestpairs <-
    function(x, y, corr_threshold=0.9, cores=1)
{
    if(n_cores(cores)==1) {
        result <- .corr_betw_matrices_unpaired_bestpairs(x, y, corr_threshold)
        result <- as.data.frame(result)
    }
    else {
        f <- function(i) {
            result <- as.data.frame(.corr_betw_matrices_unpaired_bestpairs(x[,i,drop=FALSE], y, corr_threshold))
            if(nrow(result) > 0) result$xindex <- rep(i, nrow(result))
            result
        }
        result <- cluster_lapply(cores, 1:ncol(x), f)
        result <- do.call("rbind", result)
    }

    colnames(result) <- c("corr", "xindex", "yindex")
    result$xindex <- as.integer(result$xindex)
    result$yindex <- as.integer(result$yindex)

    cbind(result,
          xcol=colnames(x)[result[,2]],
          ycol=colnames(y)[result[,3]],
          stringsAsFactors=FALSE)
}

corr_betw_matrices_unpaired_all <-
    function(x, y, cores=1)
{
    if(n_cores(cores)==1) {
        result <- .corr_betw_matrices_unpaired_all(x, y)
    }
    else {
        f <- function(i) .corr_betw_matrices_unpaired_all(x[,i,drop=FALSE], y)
        result <- cluster_lapply(cores, 1:ncol(x), f)
        result <- matrix(unlist(result), ncol=ncol(y), byrow=TRUE)
    }

    dimnames(result) <- list(colnames(x), colnames(y))
    result
}
