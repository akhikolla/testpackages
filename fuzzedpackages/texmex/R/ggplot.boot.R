ggbootdensplots <- function(x, denscol, histcol, linecol){
    v <- colnames(x$replicates)
    n <- length(v)
    p <- vector("list", length=n)

    for (i in 1:n){
        d <- data.frame(x$replicates[, i])
        names(d) <- "x"
        p[[i]] <- ggplot(data=d, aes(x=x)) +
                     geom_density(fill=denscol,colour=denscol) +
                     geom_histogram(aes(y=..density..),fill=histcol,bins=20,alpha=0.5) +
                     scale_x_continuous(v[i]) +
                     scale_y_continuous("") + 
                     geom_vline(xintercept=coef(x$map)[i], col=linecol)
    }
    p
}

#' Diagnostic plots for the replicate estimated parameter values in an evmBoot object
#' @param data An object of class 'evmBoot'.
#' @param denscol Colour for the densities. Defaults to 'light blue'.
#' @param histcol Colour for the histograms. Defaults to 'dark blue'.
#' @param linecol Colour for the point estimate lines. Decaults to 'orange'.
#' @param plot.it Whether or not to actually print the plots. Defaults
#'     to \code{plot.it=TRUE}.  If \code{plot.it=FALSE}, you might
#'     want to control the layout. Do this with
#'     \code{do.call("grid.arrange", c(plots, ncol=2))}, for example,
#'     where \code{plots} is the objected returned by
#'     \code{ggplot.evmBoot}.
#' @param mapping,environment ignored
#' @param ... Additional arguments to \code{ggplot}, currently unused.
#' @aliases ggbootdensplots
#' @keywords hplot
#' @method ggplot evmBoot
#' @export
ggplot.evmBoot <- function(data=NULL, mapping, denscol="light blue", histcol="dark blue", linecol="orange", 
                           plot.it=TRUE,
                          ..., environment){
    res <- ggbootdensplots(data, denscol=denscol, histcol=histcol, linecol=linecol)
    if (plot.it) do.call("grid.arrange", c(res, ncol=ncol(data$replicates)))
    invisible(res)
}