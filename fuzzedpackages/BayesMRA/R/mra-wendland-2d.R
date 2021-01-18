#' Code to construct the mutli-resolution sparse basis function representation for fitting spatial processes
#'
#' @param locs The location variables in 2 dimensions over which to construct the basis function representation
#' @param M The number of resolutions.
#' @param n_coarse_grid The number of basis functions in one direction (e.g. `n_coarse_grid = 10` results in a \eqn{10 \times 10}{10x10} course grid which is further extended by the number of additional padding basis functions given by `n_padding`.
#' @param n_padding The number of additional boundary points to add on each boundary. For example, n_padding = 5 will add 5 boundary knots to the both the left  and right side of the grid).
#' @param n_neighbors The expected number of neighbors for each interior basis function. This determines the basis radius parameter.
#' @param use_spam is a boolean flag to determine whether the output is a list of `spam::spam` matrix objects (`use_spam = TRUE`) or a an \eqn{n \times n}{n x n} sparse Matrix of class `Matrix::dgCMatrix` `use_spam = FALSE` (see spam and Matrix packages for details).
#'
#' @importFrom fields rdist
#' @importFrom spam spam
#' @importFrom Matrix Matrix
#' @return A list of objects including the MRA knots locations `locs_grid`,
#' the Wendland basis representation matrix `W` at the observed locations,
#' the basis radius `radius`, the numbers of resolutions `M`,
#' the number of expected neighbors in the interior of each grid `n_neighbors`,
#' the number of interior basis functions in one direction `n_coarse_grid`,
#' the number of additional padding basis functions given by `n_padding`,
#' and the setting `use_spam` which determines whether the MRA output uses the `spam` format.
#'
#' @examples
#' set.seed(111)
#' locs <- matrix(runif(20), 10, 2)
#' MRA <- mra_wendland_2d(locs, M = 2, n_coarse_grid = 4)
#' ## plot the MRA grid at different resolutions
#' layout(matrix(1:2, 1, 2))
#' plot(MRA$locs_grid[[1]])
#' plot(MRA$locs_grid[[2]])
#'
#' @export
mra_wendland_2d <- function(
    locs,
    M             = 4,
    n_coarse_grid = 10,
    n_padding     = 5L,
    n_neighbors   = 68,
    # n_max_fine_grid = 2^12,
    # radius      = 25,
    use_spam      = TRUE
) {
    ##
    ## check inputs
    ##

    if (is.null(nrow(locs))) {
        stop("locs must be a numeric matrix with N rows and 2 columns")
    }
    N <- nrow(locs)

    if (!is_numeric_matrix(locs, N, 2)) {
        stop("locs must be a numeric matrix with N rows and 2 columns")
    }
    if (!is_positive_integer(M, 1)) {
        stop("the number of resolutions M must be a positive integer")
    }
    if (!is_positive_integer(n_neighbors, 1)) {
        stop("n_neighbors must be a positive integer")
    }
    if (!is_positive_integer(n_coarse_grid, 1)) {
        stop("n_coarse_grid must be a positive integer")
    }
    if (!is_positive_integer(n_padding, 1)) {
        stop("n_padding must be a positive integer")
    }
    if (!is.logical(use_spam) || length(use_spam) != 1 || is.na(use_spam)) {
        stop("use_spam must be either TRUE or FALSE")
    }


    N <- nrow(locs)
    ## Assign as many gridpoints (approximately) as data
    # n_grid <- ceiling(sqrt(N / 2^(M:1 - 1)))

    ## finest grid is the smaller of n_max_fine_grid
    ## or the largest power of 2 larger than N
    # n_grid <- ceiling(min(n_max_fine_grid, 2^ceiling(log(N, base = 2)))^(0.5) / 2^(M:1 - 1))
    n_grid <- n_coarse_grid * 2^(1:M - 1)
    if (min(n_grid) < 4) {
        stop("There are too many resolutions to form a reliable grid. Reduce M and try again.")
    }

    ## define radius so that each basis function covers approximately 5 neighbors
    area_domain      <- diff(range(locs[, 1])) * diff(range(locs[, 2]))
    density_domain   <- max(n_grid^2) / area_domain
    radius_fine_grid <- sqrt(n_neighbors / (density_domain * base::pi))
    # radius           <- radius_fine_grid * (2^(1:M - 1))^2
    radius           <- radius_fine_grid * (2^(M:1 - 1))
    # radius <- radius / 2^(1:M - 1)

    ## generate a set of grid locations for the Wendland basis
    locs_grid    <- vector(mode = "list", length = M)
    out          <- vector(mode = "list", length = M)
    W            <- vector(mode = "list", length = M)

    for (m in 1:M) {
        ## right now assuming a 2D grid -- can generalize later
        seq_x <- seq(
            range(locs[, 1])[1],
            range(locs[, 1])[2],
            length.out = n_grid[m]
        )
        seq_y <- seq(
            range(locs[, 2])[1],
            range(locs[, 2])[2],
            length.out = n_grid[m]
        )
        delta_x <- seq_x[2] - seq_x[1]
        delta_y <- seq_y[2] - seq_y[1]
        seq_x <- c(
            min(seq_x) - delta_x * (n_padding:1),
            seq_x,
            max(seq_x) + delta_x * (1:n_padding)
        )
        seq_y <- c(
            min(seq_y) - delta_y * (n_padding:1),
            seq_y,
            max(seq_y) + delta_y * (1:n_padding)
        )

        locs_grid[[m]] <- expand.grid(seq_x, seq_y)

        D <- rdist(locs, locs_grid[[m]])
        if (use_spam) {
            ## use the spam sparse matrix package
            W[[m]] <- spam(c(wendland_basis(D, radius[m])), nrow = nrow(D), ncol = ncol(D))
        } else {
            ## use the Matrix sparse matrix package
            W[[m]] <- Matrix(wendland_basis(D, radius[m]), sparse = TRUE)
        }
    }

    out <- list(
        locs_grid     = locs_grid,
        W             = W,
        radius        = radius,
        M             = M,
        n_neighbors   = n_neighbors,
        n_coarse_grid = n_coarse_grid,
        n_padding     = n_padding,
        use_spam      = use_spam
    )

    class(out) <- "mra_wendland_2d"

    return(out)
}
