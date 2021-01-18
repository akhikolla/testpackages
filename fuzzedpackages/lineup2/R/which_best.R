#' Determine which individual has smallest distance to each individual
#'
#' For each individual represented in a distance matrix, find the
#' individual giving the smallest entry (with NAs for individuals
#' present in only the rows or only the columns).
#'
#' @param d A distance matrix
#' @param dimension Whether to get the minimum by row or by column
#' @param get_min If TRUE, get the minimum; if FALSE, get the maximum
#'
#' @return A vector with **all** distinct individuals, with the
#' character string labels for the individuals giving the minimum
#' (or maximum) value by row or column. We include all individuals
#' so that the results are aligned with the results of
#' [get_self()].
#'
#' @examples
#' # align rows in the provided dataset, lineup2ex
#' aligned <- align_matrix_rows(lineup2ex$gastroc, lineup2ex$islet)
#' # find correlated columns
#' selected_genes <- (corr_betw_matrices(aligned[[1]], aligned[[2]], "paired") > 0.75)
#' # calculate correlation between rows
#' similarity <- corr_betw_matrices(t(lineup2ex$gastroc[,selected_genes]),
#'                                  t(lineup2ex$islet[,selected_genes]), "all")
#' # which sample gives maximum value by row
#' best_byrow <- which_best(similarity, get_min=FALSE)
#'
#' # which sample gives maximum value by column
#' best_bycol <- which_best(similarity, get_min=FALSE, dimension="column")
#'
#' @seealso [get_best()], [get_self()], [get_2ndbest()], [which_2ndbest()]
#'
#' @importFrom stats setNames
#' @export
which_best <-
    function(d, dimension=c("row", "column"), get_min=TRUE)
{
    dimension <- match.arg(dimension)

    if(get_min) f <- which.min
    else f <- which.max

    rn <- rownames(d)
    cn <- colnames(d)
    if(is.null(rn) || is.null(cn))
        stop("Input matrix must have both row and column names")

    # distinct individuals
    ind <- unique(c(rn, cn))

    # pull out the best distances
    if(dimension=="row") result <- setNames(colnames(d)[apply(d, 1, f)], rownames(d))
    else result <- setNames(rownames(d)[apply(d, 2, f)], colnames(d))

    # paste into a vector with all individuals
    full_result <- setNames(rep(NA, length(ind)), ind)
    full_result[names(result)] <- result

    full_result
}
