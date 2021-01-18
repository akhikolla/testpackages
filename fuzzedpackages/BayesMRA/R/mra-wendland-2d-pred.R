#' Code to construct the mutli-resolution sparse basis function representation for fitting spatial processes
#'
#' @param locs The location variables in 2 dimensions over which to construct the basis function representation in the fitting stage.
#' @param locs_pred The location variables in 2 dimensions over which to construct the basis function representation in the prediction stage.
#' @param MRA The multi-resolution basis expansion at the observed locations. This object is the output of `mra_wendland-2d()` and is of class "mra_wendland_2d".
#' @param use_spam is a boolean flag to determine whether the output is a list of spam matrix objects (`use_spam = TRUE`) or a an \eqn{n \times n}{n x n} sparse Matrix of class "dgCMatrix" `use_spam = FALSE` (see spam and Matrix packages for details).
#'
#' @importFrom fields rdist
#' @importFrom spam spam
#' @importFrom Matrix Matrix
#' @return A list of objects including the MRA knots locations `locs_grid`,
#' the Wendland basis representation matrix `W_pred` at the prediction locations, and the basis radius `radius`
#'
#' @examples
#' set.seed(111)
#' locs <- matrix(runif(20), 10, 2)
#' locs_pred <- matrix(runif(20), 10, 2)
#' MRA <- mra_wendland_2d(locs, M = 2, n_coarse_grid = 4)
#' MRA_pred <- mra_wendland_2d_pred(locs, locs_pred, MRA)
#'
#' ## plot the MRA prediction grid at different resolutions
#' layout(matrix(1:2, 1, 2))
#' plot(MRA_pred$locs_grid[[1]])
#' plot(MRA_pred$locs_grid[[2]])
#'
#' @export
mra_wendland_2d_pred <- function(
    locs,
    locs_pred,
    MRA,
    use_spam = TRUE
) {

    N      <- nrow(locs)
    N_pred <- nrow(locs_pred)

    if (class(MRA) != "mra_wendland_2d")
        stop('MRA must be of class "mra_wendland_2d"')

    if (!is_numeric_matrix(locs, N, 2)) {
        stop("locs must be a numeric matrix with N rows and 2 columns")
    }
    if (!is_numeric_matrix(locs_pred, N_pred, 2)) {
        stop("locs_pred must be a numeric matrix with N rows and 2 columns")
    }

    if (!is.logical(use_spam) || length(use_spam) != 1 || is.na(use_spam)) {
        stop("use_spam must be either TRUE or FALSE")
    }


    M      <- length(MRA$radius)


    W_pred       <- vector(mode = "list", length = M)

    for (m in 1:M) {

        D_pred <- rdist(locs_pred, MRA$locs_grid[[m]])
        if (use_spam) {
            ## use the spam sparse matrix package
            W_pred[[m]] <- spam(c(wendland_basis(D_pred, MRA$radius[m])), nrow = nrow(D_pred), ncol = ncol(D_pred))
        } else {
            ## use the Matrix sparse matrix package
            W_pred[[m]] <- Matrix(wendland_basis(D_pred, MRA$radius[m]), sparse = TRUE)
        }
    }

    out <- list(
        locs_grid = MRA$locs_grid,
        W_pred    = W_pred,
        radius    = MRA$radius
    )

    class(out) <- "mra_wendland_2d_pred"

    return(out)
}
