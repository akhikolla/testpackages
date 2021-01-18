#' Align the rows of two matrices
#'
#' Align the rows of two matrices using their row names, omitting rows
#' that are not present in both.
#'
#' @param x A matrix
#' @param y Another matrix
#'
#' @return A list with the input `x` and `y` matrices, with the rows
#' aligned using their names. Rows not in both matrices are
#' omitted.
#'
#' @examples
#' # using the provided lineup2ex data (a list of two matrices)
#' # reduces to the common rows and puts the rows in the same order
#' # (using the row names)
#' aligned <- align_matrix_rows(lineup2ex$gastroc, lineup2ex$islet)
#'
#' @seealso [align_matrix_cols()]
#' @export
align_matrix_rows <-
    function(x, y)
{
    rnx <- rownames(x)
    rny <- rownames(y)

    if(is.null(rnx) || is.null(rny)) return(list(x=x,y=y))

    if(length(unique(rnx)) != length(rnx) ||
       length(unique(rny)) != length(rny))
        stop("The row names aren't distinct.")

    common <- rnx[rnx %in% rny]

    # return just the rows that are in common
    list(x=x[common,,drop=FALSE], y=y[common,,drop=FALSE])
}

#' Align the columns of two matrices
#'
#' Align the columns of two matrices using their column names,
#' omitting columns that are not present in both.
#'
#' @param x A matrix
#' @param y Another matrix
#'
#' @return A list with the input `x` and `y` matrices, with the
#' columns aligned using their names. Columns not in both matrices
#' are omitted.
#'
#' @examples
#' # using the provided lineup2ex data (a list of two matrices)
#' # reduces to the common columns and puts the columns in the same order
#' # (using the column names)
#' aligned <- align_matrix_cols(lineup2ex$gastroc, lineup2ex$islet)
#'
#' @seealso [align_matrix_rows()]
#' @export
align_matrix_cols <-
    function(x, y)
{
    cnx <- colnames(x)
    cny <- colnames(y)

    if(is.null(cnx) || is.null(cny)) return(list(x=x,y=y))

    if(length(unique(cnx)) != length(cnx) ||
       length(unique(cny)) != length(cny))
        stop("The column names aren't distinct.")

    common <- cnx[cnx %in% cny]

    # return just the cols that are in common
    list(x=x[,common,drop=FALSE], y=y[,common,drop=FALSE])
}
