#' Print function for S3 class "primal"
#'
#' Print the information about the model usage, the parameter path, degree of freedom of the solution path.
#'
#' @param x An object with S3 class \code{"primal"}.
#' @param \dots System reserved (No specific usage)
#' @seealso \code{\link{Dantzig_solver}}, \code{\link{SparseSVM_solver}}
#' @export
print.primal <- function(x, ...) {
    cat("\n*******Parametric Simplex Method solving ")
    cat(x$type, "problem*********\n")
    cat("iteration times = ", x$iterN, "\n")
    cat("lambda list:\n")
    print(signif(x$lambda, digits = 5))
    cat("Degree of freedom:", x$df[1], "----->", x$df[x$iterN], "\n")
    if (units.difftime(x$runtime) == "secs")
        unit <- "secs"
    if (units.difftime(x$runtime) == "mins")
        unit <- "mins"
    if (units.difftime(x$runtime) == "hours")
        unit <- "hours"
    cat("Runtime:", x$runtime, " ", unit, "\n")
}

