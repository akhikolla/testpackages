#' calculate the Wendland basis function
#'
#' @param d The distance over which to calculate the Wendland basis
#' @param radius The effective radius over which the Wendland basis is defined
#'
#' @return The output of the Wendland basis applied to the distance `d` for a given radius `radius`.
#'
#' @examples
#' layout(matrix(1:2, 1, 2))
#' curve(wendland_basis(sqrt(x^2), radius = 1), from = -2, to = 2)
#' curve(wendland_basis(sqrt(x^2), radius = 2), from = -2, to = 2)
#'
#' @export
#'
wendland_basis <- function(d, radius) {
    if (any(is.na(d))) {
        stop("d must not contain missing values")
    }
    if (any(d < 0)) {
        stop("d must be nonnegative")
    }
    if (!is_positive_numeric(radius, 1)) {
        stop("radius must be a single positive numeric value")
    }

    d_rad <- d / radius
    return(((1 - d_rad)^6 * (35 * d_rad^2 + 18 * d_rad + 3)) / 3 * (d_rad < 1))
}
