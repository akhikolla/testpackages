#' Get self-self distance
#'
#' For each individual represented in a distance matrix, pull the
#' self-self entry (with NAs for individuals present in only the rows
#' or only the columns).
#'
#' @param d A distance matrix
#'
#' @return A vector with all distinct individuals, with the self-self
#' values
#'
#' @examples
#' # align rows in the provided dataset, lineup2ex
#' aligned <- align_matrix_rows(lineup2ex$gastroc, lineup2ex$islet)
#' # find correlated columns
#' selected_genes <- (corr_betw_matrices(aligned[[1]], aligned[[2]], "paired") > 0.75)
#' # calculate correlation between rows
#' similarity <- corr_betw_matrices(t(lineup2ex$gastroc[,selected_genes]),
#'                                  t(lineup2ex$islet[,selected_genes]), "all")
#' # pull out the self-self similarities
#' self <- get_self(similarity)
#'
#' @seealso [get_best()], [get_2ndbest()], [get_nonself()]
#'
#' @export
get_self <-
    function(d)
{
    rn <- rownames(d)
    cn <- colnames(d)
    if(is.null(rn) || is.null(cn))
        stop("Input matrix must have both row and column names")

    # distinct individuals
    ind <- unique(c(rn, cn))

    # pull out the self-self distances
    vapply(ind, function(a) ifelse(a %in% rn && a %in% cn, d[a,a], NA), 1, USE.NAMES=TRUE)
}
