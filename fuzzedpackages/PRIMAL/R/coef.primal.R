#' Coef function for S3 class "primal"
#'
#' Print the estimated solution correspond to a specific parameter.
#'
#' @param object An object with S3 class \code{"primal"}.
#' @param n The index of the wanted parameter.
#' @param \dots System reserved (No specific usage)
#' @seealso \code{\link{Dantzig_solver}}, \code{\link{SparseSVM_solver}}
#' @export
coef.primal <- function(object, n = NULL, ...) {
    if (is.null(n))
        n <- object$iterN
    coef<-list(
        problem = object$type,
        index = n,
        lambda = object$lambda[n],
        df = object$df[n],
        beta = object$beta[, n])
    if(object$type == "SparseSVM"){
        coef['beta0'] = object$beta0[n]
        }
    return(coef)
    
}
