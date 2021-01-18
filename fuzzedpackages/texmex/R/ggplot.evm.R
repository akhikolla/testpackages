ggplot.qqevm <- function(data=NULL, mapping, xlab, ylab, main,
                       ylim = "auto",
                       ptcol="blue",
                       col="light blue",
                       fill="orange",
                       ..., environment){

    if (missing(xlab) || is.null(xlab)) { xlab <- "Model" }
    if (missing(ylab) || is.null(ylab)) { ylab <- "Empirical" }
    if (missing(main) || is.null(main)) { main <- "Quantile Plot" }

    limits <- if (is.null(data$sim)) NULL else range(data$sim, data$dat)

    d <- data.frame(x = data$ModPoints, y = sort(data$dat))

    poly <- data.frame(x = c(data$ModPoints,rev(data$ModPoints)),
                       y = c(data$sim[1, ], rev(data$sim[2, ])))

    p <- ggplot(poly, aes(x, y)) +
             geom_polygon(fill=fill, alpha=.5) +
             geom_abline(intercept=0, slope=1, color=col) +
             geom_point(data=d, aes(x, y), color=ptcol) +
             ggtitle(main) +
             scale_x_continuous(xlab) +
             scale_y_continuous(ylab)
    p
}

ggplot.ppevm <- function(data=NULL, mapping, xlab, ylab,  main,
                         ptcol="blue", col="light blue", fill="orange", ..., environment){

    if (missing(xlab) || is.null(xlab) ) { xlab <- "Model" }
    if (missing(ylab) || is.null(ylab) ) { ylab <- "Empirical" }
    if (missing(main) || is.null(main) ) { main <- "Probability Plot" }

    d <- data.frame(x=data$ModPoints, y=data$pfun(sort(data$dat), data$a, data$model))
    poly <- data.frame(x = c(data$ModPoints, rev(data$ModPoints)),
                       y = c(data$sim[1, ], rev(data$sim[2, ])))

    p <- ggplot(poly, aes(x, y)) +
             geom_polygon(fill=fill, alpha=.5) +
             geom_abline(intercept=0, slope=1, color=col) +
             geom_point(data=d, aes(x, y), color=ptcol) +
             ggtitle(main) +
             scale_x_continuous(xlab) +
             scale_y_continuous(ylab)

    p
}

ggplot.hist.evmOpt <- function(data, mapping, xlab=NULL, ylab=NULL, main=NULL,
                               ptcol="orange", col="dark blue", fill="light blue", ..., environment){

    if (missing(xlab) || is.null(xlab)) { xlab <- "Data" }
    if (missing(ylab) || is.null(ylab)) { ylab <- "" }
    if (missing(main) || is.null(main)) { main <- "Histogram and density" }

    d <- data.frame(x=data$dat, y=rep(0, length(data$dat)))
    dens <- data.frame(x=c(data$dens[, 1], rev(data$dens[, 1])),
                       y=c(data$dens[, 2], rep(0, nrow(data$dens))))

    p <- ggplot(dens, aes(x, y)) +
             geom_polygon(fill=fill) +
             geom_histogram(data=d, aes(x=x, y=..density..), breaks=data$breaks, fill=col, alpha=.5) +
             geom_rug(data=d, aes(x, y), sides="b", color=ptcol) +
             scale_x_continuous(xlab) +
             scale_y_continuous(ylab) +
             ggtitle(main)
    p
}

#' @export
ggplotrl <-
function(data, mapping, alpha = .050,
         xlab, ylab, main,
         ptcol = "blue", col = "light blue", fill = "orange",
         ..., environment){

    wh <- sapply(data$data$D, ncol)
    if (any(wh > 1)){
      stop("use plot method for object returned by predict.evmOpt to see rl plots if covariates in model")
    }
    if (missing(xlab) || is.null(xlab)) { xlab <- "Return period" }
    if (missing(ylab) || is.null(ylab)) { ylab <- "Return level" }
    if (missing(main) || is.null(main)) { main <- "Return Level Plot" }

    wh <- getPlotRLdata(data, alpha, NULL)

    o <- order(wh$m) # in case the return period are not in ascending order.
    wh$m <- wh$m[o]
    wh$xm <- wh$xm[o,]

    n <- wh$n; ly <- length(wh$xdat); xdat <- wh$xdat
    d <- data.frame(x= 1 / (1 - ((n - ly + 1):n) / (n + 1)),
                    y=sort(xdat))
    ln <- data.frame(x=wh$m, y=wh$xm[, 1])
    poly <- data.frame(x=c(wh$m, rev(wh$m)), y=c(wh$xm[, 2], rev(wh$xm[, 3])))

    p <- ggplot(poly, aes(x, y)) +
             geom_polygon(fill=fill, alpha=.5) +
             geom_line(data=ln, aes(x, y), color=col) +
             geom_point(data=d, aes(x, y), color=ptcol) +
             scale_x_continuous(xlab, trans="log10") +
             scale_y_continuous(ylab) +
             ggtitle(main)
    p
}

#' Diagnostic plots for an evm object
#'
#' Create and display diagnostic plots for an evm object. See \code{\link{plot.evmOpt}}
#' for further details on what is being plotted.
#' @aliases ggplot.ppevm ggplot.qqevm ggplot.hist.evmOpt
#'     ggplot.evmOpt, ggplotrl
#' @param data An object of class \code{evm}.
#' @param which Which plots to produce. Defaults to \code{which=1:4}.
#' @param main Main titles. Should have length 4.
#' @param xlab Labels for x-axes.
#' @param nsim Number of simulated datasets to create to form
#'     tolerence regions.
#' @param alpha Used to compute coverage of pointwise confidence
#'     intervals.
#' @param jitter.width Used to control the amount of horizontal
#'     jittering of points in the plots of the residuals versus
#'     covariates (when covariates are in the model).  Defaults to
#'     \code{jitter.width=0}.
#' @param span Passed to the loess smoother and defaults to
#'     \code{span=2/3}. Sometimes this choice is poor: if the loess
#'     smoother looks wrong, try \code{span=1}.
#' @param ptcol Colour for points. Defaults to \code{ptcol="blue"}.
#' @param col Colour for lines. Defaults to \code{col="light blue"}.
#' @param fill Colour for confidence regions. Defaults to
#'     \code{fill="orange"}
#' @param plot. Whether or not to display the output. Defaults to
#'     \code{plot.=TRUE}.  If the display doesn't have the desired row
#'     and column layout, the user should specify \code{plot.=FALSE},
#'     asign the output to an object, and use \code{grid.arrange} to
#'     display it.
#' @param ncol The number of columns wanted in the resulting
#'     plot. Defaults to \code{ncol=2}. This argument is passed into
#'     \code{grid.arrange}.
#' @param nrow The number of rows wanted in the resulting
#'     plot. Defaults to \code{nrow=2}. This argument is passed into
#'     \code{grid.arrange}.
#' @param ... Other arguments passed through to underlying plot
#'     functions.
#' @param mapping,environment ignored
#' @details The function attempts to arrange the plots nicely. If the
#'     output isn't what was wanted, the function returns the graphs
#'     to the user as a list so that the user can use
#'     \code{grid.arrange} directly.  Also, if you have one or more
#'     covariates in the model and the loess smoother looks wrong, try
#'     setting \code{span=1}.
#' @seealso \code{\link{plot.evmOpt}}
#' @keywords hplot
#' @method ggplot evmOpt
#' @export
ggplot.evmOpt <-
function(data, mapping, which=1:4, main=rep(NULL,4), xlab=rep(NULL,4), nsim=1000, alpha=.05, jitter.width=0,
         ptcol="blue", span=2/3, col="light blue", fill="orange", plot.=TRUE, ncol=2, nrow=2, ..., environment){
    if (!missing(main)){
        if (length(main) != 1 & length(main) != 4){
            stop("main should have length 1 or 4")
        }
        else if (length(main) == 1){ main <- rep(main, 4) }
    }

    pSymbols <- c(m="mu", f="phi", x="xi", s="sigma", k="kappa")

    pp <- qq <- rl <- h <- co <- 0

    if (all(sapply(data$data$D,ncol) == 1)){
        if (1 %in% which)
        pp <- ggplot(ppevm(data, nsim=nsim, alpha=alpha),
                     main=main[1], xlab=xlab[1],
                     ptcol=ptcol, col=col, fill=fill)
        if (2 %in% which)
        qq <- ggplot(qqevm(data, nsim=nsim, alpha=alpha),
                     main=main[2], xlab=xlab[2],
                     ptcol=ptcol, col=col, fill=fill)
        if (3 %in% which)
        rl <- ggplotrl(data, main=main[3], xlab=xlab[3], ptcol=ptcol, col=col, fill=fill)
        if (4 %in% which)
        h <- ggplot(hist.evmOpt(data), main=main[4], xlab=xlab[4])
        res <- list(pp, qq, rl, h)[which]
    } else { # Covariates in the model
        if (missing(which)){ which <- 1:3 }

        np <- length(data$data$D)
        lp <- predict(data,type="lp", unique.=FALSE)$obj[[1]]
        Which <- as.logical(apply(lp[,1:np],2,var)) # identifies which cols have covariates

        data$data$y <- resid(data)
        data$threshold <- 0
        data$coefficients <- rep(0, length(data$data$D)) # phi not sigma, so 0 not 1
        if (1 %in% which)
        pp <- ggplot(ppevm(data, nsim=nsim, alpha=alpha),
                     main=main[1], xlab=xlab[1],
                     ptcol=ptcol, col=col, fill=fill)
        if (2 %in% which)
        qq <- ggplot(qqevm(data, nsim=nsim, alpha=alpha),
                     main=main[2], xlab=xlab[2],
                     ptcol=ptcol, col=col, fill=fill)

        co <- vector("list", length=length(data$data$D[Which]))
        if (3 %in% which)
        for(i in (1:length(data$data$D))[Which]){
          ParName <- names(data$data$D[i])
          xlab <- paste("Fitted", ParName)
          d <- data.frame(lp=lp[, i], r = resid(data))

          co[[i]] <- ggplot(d, aes(lp, r)) +
                         geom_point(color=ptcol, alpha=.7, position=position_jitter(width=jitter.width)) +
                         stat_smooth(color=col, se=FALSE, span=span, method="loess") +
                         ggtitle(paste("Residuals vs fitted", ParName)) +
                         scale_x_continuous(paste("Fitted", ParName)) +
                         scale_y_continuous("Residuals", limits=range(d$r))
        }
        co <- co[!sapply(co, is.null)] # modifyList will do this

        res <- c(list(pp, qq), co)
        names(res) <- letters[1:length(res)] # stop grid.arrange getting confused

    } # Close else

    # Try to arrange the output nicely.
    #if (length(res) %% 2 == 1){
    #    blankPanel <- grid.rect(gp=gpar(col="white"))
    #    res <- c(res, list(blankPanel))
    #}

    # The loess smoother can tend to throw warnings, so suppress
    if (plot.) suppressWarnings(do.call("grid.arrange", c(res, list(ncol=ncol, nrow=nrow))))

    # Send output to the user in case it needs to be arranged
    # differently.
    invisible(res)
}
