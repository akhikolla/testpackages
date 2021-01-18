#' Estimate generalized Pareto distribution parameters over a range of values
#'
#' Estimate generalized Pareto distribution parameters over a range of values,
#' using maximum (penalized) likelihood.
#'
#' This is Stuart Coles' \code{gpd.fitrange}, as it appears in the \code{ismev}
#' package, refactored into a function that does the computations, and method
#' functions. The function uses \code{evm} internally and uses the default
#' options for that function.
#'
#' Note this function does not extend to assessing model fit when there are
#' covariates included in the model.
#'
#' @aliases gpdRangeFit print.gpdRangeFit summary.gpdRangeFit print.summary.gpdRangeFit plot.gpdRangeFit
#' ggplot.gpdRangeFit
#' @usage gpdRangeFit(data, umin=quantile(data, .05), umax=quantile(data, .95),
#' nint = 10, penalty = "gaussian", priorParameters = NULL, alpha=0.05,
#' cov="observed")
#' \method{print}{gpdRangeFit}(x, ...)
#' \method{summary}{gpdRangeFit}(object, ...)
#' \method{print}{summary.gpdRangeFit}(x, ...)
#' \method{plot}{gpdRangeFit}(x, xlab = "Threshold", ylab = NULL, main = NULL, addNexcesses=TRUE, ...)
#' \method{ggplot}{gpdRangeFit}(data, mapping, xlab="Threshold", ylab=NULL,
#' main=NULL, fill="orange", col="blue", addNexcesses = TRUE, textsize=4, ...,
#' environment)
#' @param data The data vector to be modelled.
#' @param umin The minimum threshold above which to estimate the parameters.
#' @param umax The maximum threshold above which to estimate the parameters.
#' @param nint The number of thresholds at which to perform the estimation.
#' @param penalty The type of penalty to be used in the maximum penalized
#' likelihood estimation. Should be either "gaussian" or "none". Defaults to
#' "gaussian".
#' @param priorParameters Parameters to be used for the penalty function.  See
#' the help for \code{\link{evm}} for more informaiton.
#' @param alpha 100(1 - alpha)\% confidence intervals will be plotted with the
#' point estimates. Defaults to \code{alpha = 0.05}.
#' @param cov How to compute the covariance matrix of the parameters. Defaults
#' to \code{cov = "observed"} in which case the observed information matrix is
#' used, if the \code{info} element of the \code{texmexFamily} object is
#' present. See more detailed documentation of this argument in
#' \code{\link{evm}}.
#' @param x,object Arguments to \code{print} and \code{summary} functions.
#' @param xlab Label for the x-axis.
#' @param ylab Label for the y-axis.
#' @param main The main title.
#' @param addNexcesses Annotate top axis with numbers of threshold excesses
#' arising with the corresponding values of threshold on the bottom axis.
#' @param col Colour of the line on the threshold stability plot.
#' @param fill Colour of the pointwise confidence region on the threshold
#' stability plots.
#' @param textsize Size of text on the plot (ggplot). Defaults to
#' \code{textsize=4}.
#' @param \dots Arguments to \code{plot}.
#' @param mapping,environment Not used.
#' @author Stuart Coles, Janet E Heffernan, Harry Southworth
#' @seealso \code{\link{evm}}
#' @keywords models
#' @examples
#'
#' par(mfrow=c(1,2))
#' plot(gpdRangeFit(rain))
#'
#' @export gpdRangeFit
gpdRangeFit <-
function (data, umin=quantile(data, .05), umax=quantile(data, .95),
          nint = 10,
          penalty="gaussian", priorParameters=NULL, alpha=.05,
          cov="observed") {

    m <- s <- hi <- lo <- matrix(0, nrow = nint, ncol = 2)
    u <- seq(umin, umax, length = nint)
    qz <- qnorm(1-alpha/2)
    for (i in 1:nint) {
        z <- evm(data, th=u[i], penalty=penalty, priorParameters=priorParameters, cov=cov)
        m[i, ] <- z$coefficients
        m[i, 1] <- exp(m[i, 1])

        m[i, 1] <- m[i, 1] - m[i, 2] * u[i]

        d <- matrix(c(1, -u[i]), ncol = 1)
        v <- t(d) %*% z$cov %*% d
        s[i, ] <- sqrt(diag(z$cov))
        s[i, 1] <- sqrt(v)

        hi[i, ] <- m[i, ] + qz * s[i, ]
        lo[i, ] <- m[i, ] - qz * s[i, ]
    }
    res <- list(th=u, par=m , hi=hi, lo=lo, data=data)
    oldClass(res) <- 'gpdRangeFit'
    res
}

#' @export
print.gpdRangeFit <- function(x, ...){
    sc <- cbind(threshold=x$th, phi=x$par[, 1], lo=x$lo[, 1], hi=x$hi[, 1])
    sh <- cbind(threshold=x$th, xi=x$par[, 2], lo=x$lo[, 2], hi=x$hi[, 2])

    cat("\nLog-scale parameter\n---------------\n")
    print(sc)
    cat("\nShape parameter\n---------------\n")
    print(sh)
    invisible(x)
}

#' @export
summary.gpdRangeFit <- function(object, ...){
    sc <- cbind(threshold=object$th, phi=object$par[, 1], lo=object$lo[, 1], hi=object$hi[, 1])
    sh <- cbind(threshold=object$th, xi=object$par[, 2], lo=object$lo[, 2], hi=object$hi[, 2])

    res <- list(phi=summary(sc), xi=summary(sh))
    oldClass(res) <- "summary.gpdRangeFit"
    res
}

#' @export
print.summary.gpdRangeFit <- function(x, ...){
    list(phi=x$phi, xi=x$xi)
    invisible(x)
}

#' @export
plot.gpdRangeFit <- function(x, xlab="Threshold", ylab=NULL,
                             main=NULL, addNexcesses=TRUE, ...){
    #############################################################
    ## Get axis labels and titles
    if (missing(ylab)){
        ylab <- c("scale", "shape")
    }
    else if (length( ylab ) != 2){
        stop("length of ylab should be 2")
    }

    if (!missing(main) && length(main) != 2){
        stop("length of main should be 2")
    }
    ##
    #############################################################

    data <- x$data

    for (i in 1:2) {
        yl <- range(x$hi[, i], x$lo[, i])
        plot(x$th, x$par[, i], ylim = yl, type = "b",
             xlab=xlab, ylab=ylab[i], main=main[i], ...)
        for (j in 1:length(x$th)){
            lines(c(x$th[j], x$th[j]), c(x$hi[j, i], x$lo[j, i]))
        }
        if(addNexcesses){
          axis(3,at=axTicks(1),labels=sapply(axTicks(1),function(u)sum(data>u)),cex=0.5)
          mtext("# threshold excesses")
        }
    }
    invisible()
}

