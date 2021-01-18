#' Diagnostic plots for an declustered object
#'
#' Create and display diagnostic plots for a declustered object.
#' @param data An object of class \code{declustered} or \code{extremalIndex}.
#' @param ptcol Colour for points. Defaults to \code{ptcol="blue"}.
#' @param col Colour for lines. Defaults to \code{col="light blue"}.
#' @param plot. Whether or not to display the output. Defaults to \code{plot.=TRUE}.
#' @param ... Other arguments passed through to underlying plot functions.
#' @param mapping Not used.
#' @param environment Not used.
#' @param xlab Label for the x-axis.
#' @param ylab Label for the y-axis.
#' @param main Plot title.
#' @export
ggplot.declustered <- function(data=NULL, mapping, xlab, ylab, main,
                              ptcol=c("blue","orange"),
                              col="light blue",
                              plot.=TRUE,
                              ..., environment){
    
    if (missing(xlab) || is.null(xlab)) { xlab <- "" }
    if (missing(ylab) || is.null(ylab)) { ylab <- "Data" }
    if (missing(main) || is.null(main)) { main <- ""}
    
    d <- data.frame(x = 1:length(data$y), y = data$y)

    p <- ggplot(d, aes(x, y)) +
        geom_point(data=d, aes(x, y),col=ptcol[1],...) +
        geom_hline(yintercept=data$threshold)
    
    for(i in 1:length(data$sizes)){
        d <- data.frame(x=data$exceedanceTimes[data$clusters == i], y = data$thExceedances[data$clusters == i])
        p <- p + geom_point(data=d, aes(x, y), color=ptcol[2],...)
    }
    p <- p+ggtitle(main) +
            scale_x_continuous(xlab) +
            scale_y_continuous(ylab)
    p
}

#' @rdname ggplot.declustered
#' @export
ggplot.extremalIndex <- function(data=NULL, mapping, xlab, ylab, main,
                                 ptcol="blue",
                                 col="light blue",
                                 plot.=TRUE,
                                 ..., environment){
    
    if (missing(xlab) || is.null(xlab)) { xlab <- "Standard Exponential Quantiles" }
    if (missing(ylab) || is.null(ylab)) { ylab <- "Interexceedance Times" }
    if (missing(main) || is.null(main)) { main <- paste("Threshold=",data$threshold)}
    
    NormInterExceedTimes <- data$interExceedTimes * data$thExceedanceProb
    StdExpQuantiles <- qexp(ppoints(NormInterExceedTimes))
    Theta <- data$EIintervals
    
    d <- data.frame(x = StdExpQuantiles, y = sort(NormInterExceedTimes))
    
    p <- ggplot(d, aes(x, y)) +
        geom_vline(xintercept=qexp(1-Theta)) +
        geom_abline(slope=1/Theta, intercept=-qexp(1-Theta)/Theta) +
        geom_point(data=d, aes(x, y), color=ptcol) +
        ggtitle(main) +
        scale_x_continuous(xlab) +
        scale_y_continuous(ylab)
    p
}
