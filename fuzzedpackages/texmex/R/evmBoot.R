#' Bootstrap an evmOpt fit
#'
#' This runs a parametric bootstrap simulating from an optimized
#' model.
#'
#' @param o a fit \code{evmOpt} object
#' @param R the number of parametric bootstrap samples to run
#' @param trace the frequency of trace messages
#' @param cores The number of coresto use when bootstrapping. Defaults
#'     to \code{cores=NULL} and the function guesses how many cores
#'     are available and uses them all.
#' @param export Character vector of names of variables to export. See the
#'   help file for \code{parallel::export}. Defaults to \code{export = NULL}
#'   and most users will never need to use it.
#' @param theCall (for internal use)
#' @param x an \code{\link{evmBoot}} object
#' @param col colour used to fill histogram
#' @param border the colour of the border around the bars
#' @param object a \code{\link{evmBoot}} object
#' @param ... other arguments passed to internal functions
#' @return An object of class \code{evmBoot}; a list with
#'
#' \item{call}{The call to \code{evmBoot} that produced the object.}
#' \item{replicates}{The parameter estimates from the bootstrap fits.}
#' \item{map}{The fit by by maximum penalized likelihood to the original data.}
#'
#' @aliases evmBoot summary.evmBoot plot.evmBoot coef.evmBoot print.summary.evmBoot print.evmBoot
#'
#' @usage evmBoot(o, R=1000, trace=100, cores=NULL, export=NULL, theCall)
#' \method{summary}{evmBoot}(object,...)
#' \method{plot}{evmBoot}(x,col=4,border=NULL,...)
#' \method{coef}{evmBoot}(object,...)
#' \method{print}{summary.evmBoot}(x,...)
#' \method{print}{evmBoot}(x,...)
#'
#' @note It is not expected that a user will need to call
#'     this function directly; you are directed to \code{\link{evm}}.
#' @seealso \code{\link{evm}}
#' @export
evmBoot <- function(o, R=1000, trace=100, cores=NULL, export=NULL, theCall){
    if (!inherits(o, "evmOpt")){
        stop("o must be of class 'evmOpt'")
    }

    if (missing(theCall)){ theCall <- match.call() }

    d <- o$data
    param <- texmexMakeParams(coef(o), d$D)
    rng <- o$family$rng

    getCluster <- function(n){
      wh <- try(requireNamespace("parallel"))
      if (!inherits(wh, "try-error")){
        if (is.null(n)) n <- parallel::detectCores()
        if (n == 1) { NULL }
        else parallel::makeCluster(n)
      }
      else NULL
    }
    cluster <- getCluster(cores)
    on.exit(if (!is.null(cluster)){ parallel::stopCluster(cluster) })


    bfun <- function(X){
        if (X %% trace == 0){ cat("Replicate", X, "\n") }

        s <- set.seed(seeds[[X]])

        d$y <- rng(length(d$y), param, o)

        evmFit(d, o$family, th=o$threshold, prior=o$penalty,
               priorParameters=o$priorParameters,
               start=o$coefficients,
               hessian=FALSE)$par
    }
    seeds <- as.integer(runif(R, -(2^31 - 1), 2^31))


    if (!is.null(cluster) ){
      if (!is.null(export)){
        parallel::clusterExport(cl = cluster, varlist = export)
      }
      res <- t(parallel::parSapply(cluster, X=1:R, bfun))
    }
    else {
      res <- t(sapply(1:R, bfun))
    }

    if (R > 1){
      se <- apply(res, 2, sd)
      b <- apply(res, 2, mean) - coef(o)

      if (any(abs(b/se) > .25)){
        message("Ratio of bias to standard error is high")
      }
    } # Close if (R > 1)
    res <- list(call=theCall, replicates=res, map=o)
    oldClass(res) <- "evmBoot"
    res
}

#' @export
print.evmBoot <- function(x, ...){
    print(x$call)
    means <- apply(x$replicates, 2, mean)
    medians <- apply(x$replicates, 2, median)
    sds <- rep(NA, length(means))
    if (nrow(x$replicates) > 1)
      sds <- apply(x$replicates, 2, sd)
    bias <- means - x$map$coefficients
    res <- rbind(x$map$coefficients, means, bias, sds, medians)
    rownames(res) <- c("Original", "Bootstrap mean", "Bias", "SD", "Bootstrap median")
    print(res, ...)
    if (nrow(x$replicates) > 1 & any(abs(res[3,] / res[4,]) > .25)){
        message("Ratio of bias to standard error is high")
    }
    invisible(x)
}

#' @export
coef.evmBoot <- function(object, ...){
    apply(object$replicates, 2, mean)
}

#' @export
summary.evmBoot <- function(object, ...){
    means <- apply(object$replicates, 2, mean)
    medians <- apply(object$replicates, 2, median)
    sds <- apply(object$replicates, 2, sd)
    bias <- means - coef(object$map)
    res <- rbind(coef(object$map), means, bias, sds, medians)
    rownames(res) <- c("Original", "Bootstrap mean", "Bias", "SD", "Bootstrap median")

    if (any(abs(res[3,] / res[4,]) > .25)){
        message("Ratio of bias to standard error is high")
    }

    covs <- var(object$replicates)
    res <- list(call = object$call, margins=res, covariance=covs)
    oldClass(res) <- "summary.evmBoot"
    res
}

#' @export
print.summary.evmBoot <- function(x, ...){
    print(x$call)
    print(x$margins)
    cat("\nCorrelation:\n")
    print(cov2cor(x$covariance))
    invisible(x)
}

#' @export
plot.evmBoot <- function(x, col=4, border=NULL, ...){
    pfun <- function(x, col, border, xlab,...){
        d <- density(x, n=100)
        hist(x, prob=TRUE, col=col, border=border, main="", xlab=xlab, ...)
        lines(d, lwd=2, col="grey")
        rug(x)
        invisible()
    }
    for (i in 1:ncol(x$replicates)){
        pfun(x$replicates[,i], xlab=colnames(x$replicates)[i],
             col=col, border=border)
        abline(v=coef(x$map)[i], col="cyan")
    }
    invisible()
}

