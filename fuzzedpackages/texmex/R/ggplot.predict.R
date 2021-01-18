#' Plotting function for return level estimation
#'
#' @aliases ggplot.rl.evmOpt ggplot.rl.evmSim ggplot.rl.evmBoot
#' @aliases ggplot.lp.evmOpt ggplot.lp.evmSim ggplot.lp.evmBoot
#' @param data An object of class \code{rl.evmOpt}, \code{rl.evmBoot}, \code{rl.evmSim}, \code{lp.evmOpt}, \code{lp.evmBoot} or \code{lp.evmSim}, .
#' @param ptcol Colour for points. Defaults to \code{ptcol="blue"}.
#' @param col Colour for lines. Defaults to \code{col="light blue"}.
#' @param ylim Plot limits for y-axis.
#' @param fill Colour for shading polygons.
#' @param alpha Transparency.
#' @param ... Other arguments passed through to underlying plot functions.
#' @param mapping Not used.
#' @param environment Not used.
#' @param xlab Label for the x-axis.
#' @param ylab Label for the y-axis.
#' @param main Plot title.
#' @export
ggplot.rl.evmOpt <- function(data=NULL, mapping, xlab, ylab, main,
                             ylim = "auto",
                             ptcol="blue",
                             col="light blue",
                             fill="orange",
                             alpha=0.5,
                             ..., environment){
    
    p <- plot(data,plot.=FALSE)
    
    nm <- length(names(data$obj))
    nd <- dim(data$obj[[1]])[2]
    ncov <- length(unlist(data$obj)) / (nm * nd)
    if(ncov>1){
        covnames <- dimnames(p$xm)[[2]][-(1:3)]
    } else {
        covnames <- ""
    }
    
    if (missing(main) || is.null(main)) {
        main <- "Return Level Plot"
        SetMain <- TRUE
    } else {
        if(length(main) == 1){
            SetMain <- TRUE
        } else {
            SetMain <- FALSE
        }
    }
    
    Plots <- list(NULL)
    
    for(i in 1:ncov){
        xm <- t(p$xm[i,,])
        
        d <- data.frame(m=p$m, rl=xm[,1],Lower= xm[,2],Upper=xm[,3])
        cov <- xm[1,-(1:3)]

        if(SetMain){
            if(length(covnames) == 1 && covnames != ""){
                Main <- paste(main,"\n", paste(covnames,"=",signif(cov,4),collapse=", "))
            } else {
                Main <- main
            }
        } else {
            Main <- main[i]
        }

        Plots[[i]] <- ggplot(d,aes(x=m,y=rl))+
                             labs(title = Main)+
                             geom_line() +
                             scale_x_continuous(trans="log",breaks=function(x) signif(exp(seq(from=log(x[1]),to=log(x[2]),len=5)),1))+
                             geom_polygon(data= data.frame(x=c(d$m,rev(d$m)),y=c(d$Upper,rev(d$Lower))),
                                          aes(x=x,y=y),fill=fill,alpha=alpha)
            
    }
    
    if (ncov < 5) {
        res <- c(Plots); names(res) <- letters[1:length(res)] # stop grid.arrange getting confused
        do.call("grid.arrange", c(res, ncol=length(Plots)))
    } else {
        message("ggplot.rl.evmOpt produced more than 4 plots; returning invisibly.\n\nUse gridExtra::grid.arrange for plot layout.\n\n")
    }
    invisible(Plots)
}

#' @export
ggplot.rl.evmSim <- ggplot.rl.evmOpt

#' @export
ggplot.rl.evmBoot <- ggplot.rl.evmOpt

#' @export
ggplot.lp.evmOpt <- function(data=NULL, mapping, xlab, ylab, main,
                             ylim = "auto",
                             ptcol="blue",
                             col="light blue",
                             fill="orange",
                             alpha=0.5,
                             ..., environment){
    
    p <- plot(data,plot.=FALSE)

    Plots <- list(NULL)
    
    for(i in 1:length(p)){
        d <- data.frame(x=p[[i]]$x,y=p[[i]]$y[,1],Lower = p[[i]]$y[,2],Upper=p[[i]]$y[,3])

        Plots[[i]] <- ggplot(d,aes(x=x,y=y))+
            labs(x=p[[i]]$CovName,y=p[[i]]$ParName)+
            geom_line() +
            geom_polygon(data= data.frame(x=c(d$x,rev(d$x)),y=c(d$Upper,rev(d$Lower))),
                         aes(x=x,y=y),fill=fill,alpha=alpha)
        
    }
    if (length(p) < 5) {
        res <- c(Plots); names(res) <- letters[1:length(res)] # stop grid.arrange getting confused
        do.call("grid.arrange", c(res, ncol=length(Plots)))
    } else {
        message("ggplot.lp.evmOpt produced more than 4 plots; returning invisibly.\n\nUse gridExtra::grid.arrange for plot layout.\n\n")
    }
    invisible(Plots)
}

#' @export
ggplot.lp.evmSim <- ggplot.lp.evmOpt

#' @export
ggplot.lp.evmBoot <- ggplot.lp.evmOpt
