#' print a LassoGEE object
#'
#' Print a summary of the results of a LassoGEE model.
#'
#' A summary of the cross-validated fit is produced. \code{print.cv.LassoGEE(object)}
#' will print the summary includes Working Correlation and Returned Error Value.
#'
#' @aliases print.LassoGEE
#' @param x fitted 'LassoGEE' object
#' @param digits significant digits in printout
#' @param \dots additional print arguments
#' @seealso \code{LassoGEE}, and \code{cv.LassoGEE} methods.
#' @references Li, Y., Gao, X., and Xu, W. (2020). Statistical consistency for
#' generalized estimating equation with \eqn{L_1} regularization.
#' @keywords models regression
#' @method print LassoGEE
#' @export

print.LassoGEE <- function(x, digits = NULL, ...)
{
    if(is.null(digits)) digits <- options()$digits else options(digits =
                                                                digits)
    cat("\n", x$title)
    cat("\n", x$version, "\n")
    cat("\nModel:\n")
    cat(" Link:                     ", x$model$link, "\n")
    cat(" Variance to Mean Relation:",x$model$varfun,"\n")
    if(!is.null(x$model$M))
        cat(" Correlation Structure:    ", x$model$corstr, ", M =", x$
            model$M, "\n")
    else cat(" Correlation Structure:    ", x$model$corstr, "\n")
    cat("\nCall:\n")
    dput(x$call)                        #       cat("\nTerms:\n")
###        ys <- matrix(rep(as.matrix(x$id, ncol = 1), 5), ncol = 5)
    ys <- matrix(rep(matrix(x$id, ncol = 1), 5), ncol = 5)
    ys[, 2] <- x$y
    ys[, 3] <- x$linear.predictors
    ys[, 4] <- x$fitted.values
    ys[, 5] <- x$residuals
    dimnames(ys) <- list(1:length(x$y), c("ID", "Y", "LP", "fitted",
                                          "Residual")) #       cat("\nFitted Values:\n")
    cat("\nNumber of observations : ", x$nobs, "\n")
    cat("\nMaximum cluster size   : ", x$max.id, "\n")
    nas <- x$nas
    if(any(nas))
        cat("\n\nCoefficients: (", sum(nas),
            " not defined because of singularities)\n", sep = "")
    else cat("\n\nCoefficients:\n")
    print(x$betaest, digits = digits)
    cat("\nEstimated Scale Parameter: ", format(round(x$scale, digits)))
    cat("\nLambda value: ", format(round(x$lambda.value, digits)))
    cat("\nNumber of Iterations: ", x$outer.iter)
    cat("\n\nWorking Correlation[1:4,1:4]\n")
    print(x$working.correlation[1:4, 1:4], digits = digits)
    cat("\n\nReturned Error Value:\n")
    print(x$error)
    invisible(x)
}
