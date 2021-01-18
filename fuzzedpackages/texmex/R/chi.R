#' Measures of extremal dependence
#' 
#' Compute measures of extremal dependence for 2 variables.
#' 
#' Computes the functions chi and chi-bar described by Coles, Heffernan and
#' Tawn (1999). The limiting values of these functions as the quantile
#' approaches 1 give an empirical measure of the type and strength of tail
#' dependendce exhibited by the data.
#' 
#' A limiting value of ChiBar equal to 1 indicates Asymptotic Dependence, in
#' which case the limiting value of Chi gives a measure of the strength of
#' dependence in this class.  A limiting value of ChiBar of less than 1
#' indicates Asymptotic Independence in which case Chi is irrelevant and the
#' limiting value of ChiBar gives a measure of the strength of dependence.
#' 
#' The plot and ggplot methods show the ChiBar and Chi functions.  In the case of the
#' confidence interval for ChiBar excluding the value 1 for all of the largest
#' quantiles, the plot of the Chi function is shown in grey.
#' 
#' @aliases chi summary.chi plot.chi print.chi ggplot.chi print.summary.chi
#' 
#' @usage chi(data, nq = 100, qlim = NULL, alpha = 0.05, trunc = TRUE)
#' 
#' \method{summary}{chi}(object, ...)
#' 
#' \method{print}{summary.chi}(x, digits=3, ...)
#' 
#' \method{print}{chi}(x, ...)
#'
#' \method{plot}{chi}(x, show=c("Chi"=TRUE,"ChiBar"=TRUE), lty=1,
#' cilty=2, col=1, spcases=TRUE, cicol=1, xlim=c(0, 1), ylimChi =
#' c(-1, 1), ylimChiBar = c(-1, 1), mainChi = "Chi", mainChiBar =
#' "Chi Bar", xlab = "Quantile", ylabChi =
#' expression(chi(u)), ylabChiBar = expression(bar(chi)(u)),
#' ask, ...)
#' 
#' \method{ggplot}{chi}(data=NULL, mapping, xlab = "Quantile", 
#' ylab=c("ChiBar" = expression(bar(chi)(u)), "Chi" = expression(chi(u))),
#' main=c("ChiBar" = "Chi Bar",       "Chi" = "Chi"),
#' xlim = c(0, 1), ylim =list("Chi" = c(-1, 1),"ChiBar" = c(-1, 1)),
#' ptcol="blue",fill="orange",show=c("ChiBar"=TRUE,"Chi"=TRUE),
#' spcases = TRUE,plot., ..., environment)
#' 
#' @param data A matrix containing 2 numeric columns.
#' @param nq The number of quantiles at which to evaluate the dependence
#' measures.
#' @param qlim The minimum and maximum quantiles at which to do the evaluation.
#' @param alpha The size of the confidence interval to be used. Defaults to
#' \code{alpha = 0.05}.
#' @param trunc Logical flag indicating whether the estimates should be
#' truncated at their theoretical bounds.  Defaults to \code{trunc = TRUE}.
#' @param x,object An object of class \code{chi}.
#' @param digits Number of digits for printing.
#' @param show Logical, of length 2, names "Chi" and "ChiBar".  Defaults to\cr
#' \code{c("Chi" = TRUE, "ChiBar" = TRUE)}.
#' @param lty,cilty,col,cicol Line types and colours for the the estimated
#' quantities and their confidence intervals.
#' @param xlim,ylimChi,ylimChiBar Limits for the axes.
#' @param mainChi,mainChiBar Main titles for the plots.
#' @param xlab,ylabChi,ylabChiBar Axis labels for the plots.
#' @param mapping,ylab,main,ylim,ptcol,fill,environment Arguments to ggplot methods.
#' @param spcases Whether or not to plot special cases of perfect (positive and
#' negative) dependence and indpenendence. Defaults to \code{FALSE}.
#' @param plot. whether to plot to active graphics device.
#' @param ask Whether or not to ask before reusing the graphics device.
#' @param ... Further arguments to be passed to methods.
#' @return An object of class \code{chi} containing the following.
#' 
#' \item{chi}{Values of chi and their estimated upper and lower confidence
#' limits.} \item{chibar }{Values of chibar and their estimated upper and lower
#' confidence limits.} \item{quantile}{The quantiles at which chi and chi-bar
#' were evaluated.} \item{chiulb, chibarulb}{Upper and lower bounds for chi and
#' chi-bar.}
#' @note When the data contain ties, the values of chi and chibar are
#' calculated by assigning distinct ranks to tied values using the \code{rank}
#' function with argument \code{ties.method = "first"}.  This results in the
#' values of chi and chibar being sensitive to the order in which the tied
#' values appear in the data.
#' 
#' The code is a fairly simple reorganization of code written by Janet E.
#' Heffernan and Alec Stephenson and which appears in the \code{chiplot}
#' function in the \code{evd} package.
#' @author Janet E. Heffernan, Alec Stephenson, Harry Southworth
#' @seealso \code{\link{MCS}}, \code{\link{rank}}
#' @references S. Coles, J. E. Heffernan and J. A. Tawn, Dependence measures
#' for extreme values analyses, Extremes, 2, 339 -- 365, 1999.
#' 
#' A. G. Stephenson. evd: Extreme Value Distributions, R News, 2, 2002.
#' @keywords multivariate
#' @examples
#' 
#' 
#' D <- liver[liver$dose == "D",]
#' chiD <- chi(D[, 5:6])
#' par(mfrow=c(1,2))
#' ggplot(chiD)
#' 
#' A <- liver[liver$dose == "A",]
#' chiA <- chi(A[, 5:6])
#' # here the limiting value of chi bar(u) lies away from one so the chi plot is
#' # not relevant and is plotted in grey
#' ggplot(chiA) 
#' 
#' 
#' 
#' @export chi
chi <- 
  # Much of the code in here was written by Alec Stephenson
  # and comes from his 'chiplot' function in his 'evd' package.
  # Minor differences between the evd implementation and this are: 
  # use of ties.method="first" here as oppsed to the evd method 
  # which used ties.method="average"; lower bound on chibar wrong 
  # in evd package.
function (data, nq = 100, qlim = NULL, alpha = 0.05, trunc = TRUE) {
    
    theCall <- match.call()

    eps <- .Machine$double.eps^0.5

    data <- na.omit(data)
    n <- nrow(data)

    # Get the EDFs
    t.method <- "first"

    data <- cbind(rank(data[, 1],ties.method = t.method)/(n + 1), 
                  rank(data[, 2],ties.method = t.method)/(n + 1))
	
    rowmax <- apply(data, 1, max)
    rowmin <- apply(data, 1, min)

    qlim2 <- c(min(rowmax) + eps, max(rowmin) - eps)

    if (!is.null(qlim)) {
        if (qlim[1] < qlim2[1]){ stop("lower quantile limit is too low") }
        if (qlim[2] > qlim2[2]){ stop("upper quantile limit is too high") }
        if (qlim[1] > qlim[2]){ stop("lower quantile limit is less than upper quantile limit") }
    }
    else{
        qlim <- qlim2
    }

    u <- seq(qlim[1], qlim[2], length = nq)

    # HS. Replaced 2 for loops with calls to sapply

    cu <- vapply(seq_len(nq),
                 function(i, x, y){ mean(y < x[i]) },
                 0, y=rowmax, x=u )
    cbaru <- vapply(seq_len(nq), function(i, x, y){ mean(y > x[i]) },
                    0, y=rowmin, x=u )

    # Get \chi and \bar\chi
    chiu <- 2 - log(cu)/log(u) # 3.2 of Coles, Heffernan, Tawn
    chibaru <- (2 * log(1 - u))/log(cbaru) - 1 # Page 348 of Coles, Heffernan, Tawn

    # Get confidence limits
    varchi <- ((1/log(u)^2 * 1)/cu^2 * cu * (1 - cu))/n
    varchi <- qnorm(1 - alpha/2) * sqrt(varchi)

    varchibar <- (((4 * log(1 - u)^2)/(log(cbaru)^4 * cbaru^2)) * cbaru * (1 - cbaru))/n
    varchibar <- qnorm(1 - alpha/2) * sqrt(varchibar)

    chiu <- cbind(chilow = chiu - varchi,           # Lower
                  chi = chiu,                       # Point est.
                  chiupp = chiu + varchi)           # Upper

    chibaru <- cbind(chiblow = chibaru - varchibar, # Lower
                     chib = chibaru,                # Point est.
                     chibupp = chibaru + varchibar) # Upper

    chiulb <- 2 - log(pmax(2 * u - 1, 0))/log(u)

    if (trunc) {
        chiu[chiu > 1] <- 1
        chiu <- apply(chiu, 2, function(x, z){ pmax(x, z) }, z = chiulb)
        chibaru[chibaru > 1] <- 1
        chibaru[chibaru < -1] <- -1
    }

    res <- list(chi=chiu, chibar = chibaru, quantile=u, call=theCall, 
                qlim=qlim, chiulb = chiulb)
    oldClass(res) <- "chi"
    res
} # Close chi <- function


#' @export
print.chi <- function(x, ...){
    print(x$call,...)
    cat("Values of chi and chi-bar obtained and",
         length(x$quantile), "quantiles.\n")
    invisible(x)
}

#' @export
summary.chi <- function(object, ...){
    wh <- quantile(object$quantile, prob=c(.05, .5, .95))
    wh <- sapply(wh, function(i, u){
                        d <- abs(u - i)
                        min(u[d == min(d)])
                     }, u=object$quantile)

    chiQ <- object$chi[object$quantile %in% wh, 2]
    chibarQ <- object$chibar[object$quantile %in% wh, 2]

    out <- rbind(wh, chiQ, chibarQ)
    dimnames(out) <- list(c("quantile", "chi", "chi-bar"), rep("", 3))
    res <- list(out=out,call=object$call,quantile=object$quantile)
    oldClass(res) <- "summary.chi"
    res
}

#' @export
print.summary.chi <- function(x,digits=3,...){
    cat("Call: ")
    print(x$call)
    cat("Values of chi and chi-bar obtained at",
        length(x$quantile), "quantiles.\n")
    
    print(x$out, digits=digits, ...)
    invisible(x)
}

#' @export
plot.chi <- function(x, show=c("Chi"=TRUE,"ChiBar"=TRUE), lty = 1, cilty = 2, col = 1, spcases = TRUE, cicol = 1,
                     xlim = c(0, 1), ylimChi = c(-1, 1), ylimChiBar = c(-1, 1),
                     mainChi = "Chi", mainChiBar = "Chi Bar",
                     xlab = "Quantile", 
                     ylabChi = expression(chi(u)),#"Chi",
                     ylabChiBar = expression(bar(chi)(u)), #"Chi Bar",
                     ask, ...){

  lty <- c(cilty, lty, cilty)
  col <- c(cicol, col, cicol)
  nb.fig <- prod(par("mfcol"))

    if (missing(ask)){
        ask <- nb.fig < sum(show) && dev.interactive()
    }
    
  if (ask) {
     op <- par(ask = TRUE)
     on.exit(par(op))
  }
  ChiBarAsympIndep <- prod(tail(x$chibar[,3]) < 1)
  if (show["ChiBar"]) {
    matplot(x$quantile, x$chibar, type = "l", lty = lty, col = col, 
            xlim = xlim, ylim = ylimChiBar, main = mainChiBar, xlab = xlab, 
            ylab = ylabChiBar, ...)
    if (spcases) {
      segments(x$qlim[1], 0, x$qlim[2], 0, lty = 5, col = "grey")
      segments(x$qlim[1], 1, x$qlim[2], 1, lty = 5, col = "grey")
      segments(x$qlim[1],-1, x$qlim[2],-1, lty = 5, col = "grey")
    }
  }
  if (show["Chi"]) {
    if (ChiBarAsympIndep) {
      col <- "grey"
      cols <- "grey"
    } else {
      cols <- "black"
    }
  
    matplot(x$quantile, x$chi, type = "l", lty = lty, col = col, xlim = xlim, 
            ylim = ylimChi, main = mainChi, xlab = xlab, ylab = ylabChi, 
            col.lab=cols,col.main=cols,col.sub=cols,axes=FALSE,
             ...)
     box(col=cols)
     axis(1,col=cols,col.axis = cols)
     axis(2,col=cols,col.axis = cols)
             
     if (spcases) {
       segments(x$qlim[1], 0, x$qlim[2], 0, lty = 5, col = "grey")
       segments(x$qlim[1], 1, x$qlim[2], 1, lty = 5, col = "grey")
       lines(x$quantile, x$chiulb, lty = 5, col = "grey")
     }
  }

  invisible()
}

#' @export
ggplot.chi <- function(data=NULL, mapping, 
                       xlab = "Quantile", 
                       ylab=c("ChiBar" = expression(bar(chi)(u)),
                              "Chi" = expression(chi(u))),
                       main=c("ChiBar" = "Chi Bar",
                              "Chi" = "Chi"),
                       xlim = c(0, 1), 
                       ylim =list("Chi" = c(-1, 1), 
                               "ChiBar" = c(-1, 1)),
                       ptcol="blue",
                       fill="orange",
                       show=c("ChiBar"=TRUE,
                              "Chi"=TRUE),
                       spcases = TRUE,
                       plot.=TRUE,
                       ..., environment){
    
    ChiBarAsympIndep <- prod(tail(data$chibar[,3]) < 1)
    
    dat <- data.frame(qu = data$quantile, 
                      chi = data$chi[,"chi"], 
                      chibar = data$chibar[,"chib"])
    
    poly <- data.frame(qu = c(data$quantile,rev(data$quantile)),
                       chi = c(data$chi[,"chilow"],rev(data$chi[,"chiupp"])),
                       chibar = c(data$chibar[,"chiblow"],rev(data$chibar[,"chibupp"])))
    
    if (show["ChiBar"]) {
        p1 <- ggplot(dat,aes(qu,chibar)) + 
            geom_line(colour=ptcol) + 
            geom_polygon(data=poly,mapping=aes(qu,chibar),fill=fill,alpha=0.5) +
            coord_cartesian(xlim=xlim,ylim=ylim$ChiBar) +
            labs(x=xlab,y=ylab["ChiBar"],title=main["ChiBar"])
        
        if (spcases) {
            p1 <- p1 + 
                geom_line(data=data.frame(x=data$qlim,y=c(0,0)),aes(x,y),col="grey") +
                geom_line(data=data.frame(x=data$qlim,y=c(1,1)),aes(x,y),col="grey") + 
                geom_line(data=data.frame(x=data$qlim,y=c(-1,-1)),aes(x,y),col="grey") 
        }
    }
    if (show["Chi"]) {
        if (ChiBarAsympIndep) {
            ptcol <- "dark grey"
            fill <- "grey"
        }
        p2 <- ggplot(dat,aes(qu,chi)) + 
            geom_line(colour=ptcol) + 
            geom_polygon(data=poly,mapping=aes(qu,chi),fill=fill,alpha=0.5) +
            coord_cartesian(xlim=xlim,ylim=ylim$Chi) +
            labs(x=xlab,y=ylab["Chi"],title=main["Chi"])
        
        if (spcases) {
            p2 <- p2 + 
                geom_line(data=data.frame(x=data$qlim,y=c(0,0)),aes(x,y),col="grey") +
                geom_line(data=data.frame(x=data$qlim,y=c(1,1)),aes(x,y),col="grey") + 
                geom_line(data=data.frame(x=data$quantile,y=data$chiulb),aes(x,y),col="grey")
        }
    }
    
    if(plot.)gridExtra::grid.arrange(p1,p2,ncol=2)
    invisible(list(p1,p2))
}
