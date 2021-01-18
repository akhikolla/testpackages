#' Determine which individual has 2nd-smallest distance to each individual
#'
#' For each individual represented in a distance matrix, find the
#' individual giving the 2nd-smallest entry (with NAs for individuals
#' present in only the rows or only the columns).
#'
#' @param d A distance matrix
#' @param dimension Whether to get the 2nd-best by row or by column
#' @param get_min If TRUE, get the 2nd-minimum; if FALSE, get the 2nd-maximum
#'
#' @return A vector with **all** distinct individuals, with the
#' character string labels for the individuals giving the
#' 2nd-smallest (or largest) value by row or column. We include
#' all individuals so that the results are aligned with the
#' results of [get_self()].
#'
#' @examples
#' # align rows in the provided dataset, lineup2ex
#' aligned <- align_matrix_rows(lineup2ex$gastroc, lineup2ex$islet)
#' # find correlated columns
#' selected_genes <- (corr_betw_matrices(aligned[[1]], aligned[[2]], "paired") > 0.75)
#' # calculate correlation between rows
#' similarity <- corr_betw_matrices(t(lineup2ex$gastroc[,selected_genes]),
#'                                  t(lineup2ex$islet[,selected_genes]), "all")
#' # which sample gives second-biggest value by row
#' secbest_byrow <- which_2ndbest(similarity, get_min=FALSE)
#'
#' # which sample gives second-biggest value by column
#' secbest_bycol <- which_2ndbest(similarity, get_min=FALSE, dimension="column")
#'
#' @seealso [get_2ndbest()], [get_self()], [get_best()], [which_best()]
#'
#' @importFrom stats setNames
#' @export
which_2ndbest <-
    function(d, dimension=c("row", "column"), get_min=TRUE)
{
    dimension <- match.arg(dimension)

    # function that does the key work
    if(get_min) { orderf <- function(x) order(x, na.last=TRUE, decreasing=FALSE) }
    else { orderf <- function(x) order(x, na.last=TRUE, decreasing=TRUE) }
    f <- function(x, nam) nam[orderf(x)[2]]

    rn <- rownames(d)
    cn <- colnames(d)
    if(is.null(rn) || is.null(cn))
        stop("Input matrix must have both row and column names")

    # distinct individuals
    ind <- unique(c(rn, cn))

    # pull out the best distances
    if(dimension=="row") result <- setNames(apply(d, 1, f, colnames(d)), rownames(d))
    else result <- setNames(apply(d, 2, f, rownames(d)), colnames(d))

    # paste into a vector with all individuals
    full_result <- setNames(rep(NA, length(ind)), ind)
    full_result[names(result)] <- result

    full_result
}
