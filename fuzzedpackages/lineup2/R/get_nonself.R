#' Get self-nonself distances
#'
#' Return the distance matrix with all self-self distances replaced
#' with NAs (and so just containing the self-self distances).
#'
#' @param d A distance matrix
#'
#' @return The input distance matrix with all self-self distances
#' replaced with NAs.
#'
#' @examples
#' # align rows in the provided dataset, lineup2ex
#' aligned <- align_matrix_rows(lineup2ex$gastroc, lineup2ex$islet)
#' # find correlated columns
#' selected_genes <- (corr_betw_matrices(aligned[[1]], aligned[[2]], "paired") > 0.75)
#' # calculate correlation between rows
#' similarity <- corr_betw_matrices(t(lineup2ex$gastroc[,selected_genes]),
#'                                  t(lineup2ex$islet[,selected_genes]), "all")
#' # pull out the non-self similarities
#' nonself <- get_nonself(similarity)
#'
#' @seealso [get_self()], [get_best()], [get_2ndbest()]
#'
#' @export
get_nonself <-
    function(d)
{
    rn <- rownames(d)
    cn <- colnames(d)
    if(is.null(rn) || is.null(cn))
        stop("Input matrix must have both row and column names")

    # distinct individuals
    common <- rn[rn %in% cn]

    # omit self-self distances
    for(ind in common) d[ind,ind] <- NA

    d
}
