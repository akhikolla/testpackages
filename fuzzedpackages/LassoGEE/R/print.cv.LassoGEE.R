#' print a cross-validated LassoGEE object
#'
#' Print a summary of the results of cross-validation for a LassoGEE model.
#'
#' A summary of the cross-validated fit is produced. \code{print.cv.LassoGEE(object)}
#' will print the summary for a sequence of \code{lambda}.
#'
#' @aliases print.cv.LassoGEE
#' @param x fitted 'cv.LassoGEE' object
#' @param digits significant digits in printout
#' @param \dots additional print arguments
#' @seealso \code{LassoGEE}, and \code{cv.LassoGEE} methods.
#' @references Li, Y., Gao, X., and Xu, W. (2020). Statistical consistency for
#' generalized estimating equation with \eqn{L_1} regularization.
#' @keywords models regression
#' @method print cv.LassoGEE
#' @export

print.cv.LassoGEE <- function(x, digits = NULL, ...)
{
    if(is.null(digits)) digits <- options()$digits else options(digits =
                                                                digits)

    cat("\nCall:\n")
    dput(x$call)                        #       cat("\nTerms:\n")

    cat(sprintf("\n%d-fold CV results:\n", x$fold))
    print(cbind("lambda"=x$lam.vect, "Cv"=x$cv.vect))
    cat("\nOptimal tuning parameter:\n")
    optimalTuning <- c("Best lambda"=x$lam.opt)
    print(optimalTuning)

    # return object invisibly
    invisible(x)
}
