
#' Plot method for cv.vennLasso fitted objects
#'
#' @param sign.lambda Either plot against log(lambda) (default) or its negative if \code{sign.lambda = -1}.
#' @rdname plot
#' @method plot cv.vennLasso
#' @importFrom graphics abline
#' @importFrom graphics axis
#' @importFrom graphics points
#' @importFrom graphics segments
#' @importFrom graphics title
#' @export
#' @examples
#' set.seed(123)
#' 
#' dat.sim <- genHierSparseData(ncats = 3, nvars = 25,
#'                              nobs = 100, 
#'                              hier.sparsity.param = 0.5,
#'                              prop.zero.vars = 0.5,
#'                              effect.size.max = 0.25,
#'                              family = "gaussian")
#'
#' x        <- dat.sim$x
#' x.test   <- dat.sim$x.test
#' y        <- dat.sim$y
#' y.test   <- dat.sim$y.test
#' grp      <- dat.sim$group.ind
#' grp.test <- dat.sim$group.ind.test
#'
#' fit.adapt <- cv.vennLasso(x, y,
#'                           grp,
#'                           adaptive.lasso = TRUE,
#'                           nlambda        = 25,
#'                           nfolds         = 4)
#'                                      
#' plot(fit.adapt) 
#' 
plot.cv.vennLasso <- function(x, sign.lambda = 1, ...) 
{
    # compute total number of selected variables for each
    # tuning parameter 
    nzero <- apply(x$vennLasso.fit$beta[,-1,], 3, function(bb) sum(bb != 0))
    
    xlab <- expression(log(lambda))
    if(sign.lambda<0)xlab <- paste("-", xlab, sep="")
    plot.args = list(x    = sign.lambda * log(x$lambda),
                     y    = x$cvm,
                     ylim = range(x$cvup, x$cvlo),
                     xlab = xlab,
                     ylab = x$name,
                     type = "n")
    new.args <- list(...)
    if(length(new.args)) plot.args[names(new.args)] <- new.args
    do.call("plot", plot.args)
    error.bars(sign.lambda * log(x$lambda), 
               x$cvup, 
               x$cvlo, width = 0.005)
    points(sign.lambda*log(x$lambda), x$cvm, pch=20, col="dodgerblue")
    axis(side   = 3,
         at     = sign.lambda * log(x$lambda),
         labels = paste(nzero), tick=FALSE, line=0, ...)
    abline(v = sign.lambda * log(x$lambda.min), lty=2, lwd = 2, col = "firebrick3")
    abline(v = sign.lambda * log(x$lambda.1se), lty=2, lwd = 2, col = "firebrick1")
    title(x$name, line = 2.5, ...)
    
}
