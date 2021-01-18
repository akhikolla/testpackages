#' Estimate dependence parameters in a conditional multivariate extreme values
#' model over a range of thresholds.
#'
#' Diagnostic tool to aid the choice of threshold to be used for the estimation
#' of the dependence parameters in the conditional multivariate extreme values
#' model of Heffernan and Tawn, 2004.
#'
#' Dependence model parameters are estimated using a range of threshold values.
#' The sampling variability of these estimates is characterised using the
#' bootstrap.  Point estimates and bootstrap estimates are finally plotted over
#' the range of thresholds.  Choice of threshold should be made such that the
#' point estimates at the chosen threshold and beyond are constant, up to
#' sampling variation.
#'
#' @usage mexRangeFit(x, which, quantiles = seq(0.5, 0.9, length = 9),
#' start=c(.01, .01), R = 10, nPass=3, trace=10, margins = "laplace", constrain
#' = TRUE, v = 10, referenceMargin=NULL)
#' @param x An object of class \code{\link{mex}} or \code{\link{migpd}}.
#' @param which The variable on which to condition.
#' @param quantiles A numeric vector specifying the quantiles of the marginal
#' distribution of the conditioning variable at which to fit the dependence
#' model.
#' @param start See documentation for this argument in
#' \code{\link{mexDependence}}.
#' @param R The number of bootstrap runs to perform at each threshold. Defaults
#' to \code{R}=10.
#' @param nPass Argument passed to function \code{\link{bootmex}}.
#' @param trace Argument passed to function \code{\link{bootmex}}.
#' @param margins Argument passed to function \code{\link{mexDependence}}.
#' @param constrain Argument passed to function \code{\link{mexDependence}}.
#' @param v Argument passed to function \code{\link{mexDependence}}.
#' @param referenceMargin Optional set of reference marginal distributions to use for marginal transformation if the data's own marginal distribution is not appropriate (for instance if only data for which one variable is large is available, the marginal distributions of the other variables will not be represented by the available data).  This object can be created from a combination of datasets and fitted GPDs using the function \code{makeReferenceMarginalDistribution}.
#' @return NULL.
#' @author Harry Southworth, Janet E. Heffernan
#' @seealso \code{\link{mexDependence}}, \code{\link{bootmex}}
#' @references J. E. Heffernan and J. A. Tawn, A conditional approach for
#' multivariate extreme values, Journal of the Royal Statistical society B, 66,
#' 497 -- 546, 2004
#' @keywords models multivariate
#' @examples
#'
#' \donttest{
#'   w <- migpd(winter, mqu=.7)
#'   w
#'   par(mfrow=c(4,2))
#'   plot(mexRangeFit(w, which=1),main="Winter data, Heffernan and Tawn 2004",cex=0.5)
#' }
#'
#' @export mexRangeFit
mexRangeFit <-
function (x, which, quantiles=seq(0.5,0.9,length=9), start=c(.01, .01), R=10, nPass=3, trace=10,
          margins="laplace", constrain=TRUE, v=10, referenceMargin=NULL){
  if (inherits(x, "mex")){
    if( (!missing(margins))){
      warning("margins given, but already specified in 'mex' object.  Using 'mex' value")
    }
    if( (!missing(constrain))){
      warning("constrain given, but already specified in 'mex' object.  Using 'mex' value")
    }
    if( (!missing(v))){
      warning("v given, but already specified in 'mex' object.  Using 'mex' value")
    }
    if( (!missing(which))){
      warning("which given, but already specified in 'mex' object.  Using 'mex' value")
    }
    if( (!missing(referenceMargin))){
      warning("referenceMargin given, but already specified in 'mex' object.  Using 'mex' value")
    }
    constrain <- x$dependence$constrain
    v <- x$dependence$v
    which <- x$dependence$which
    margins <- x$dependence$margins
    referenceMargin <- x$referenceMargin
    x <- x[[1]]
  } else {
    if (!inherits(x, "migpd")){
      stop("object should have class mex or migpd")
    }
    if (missing(which)) {
      which <- 1
      message(paste("Missing 'which'. Conditioning on", names(x$models)[which], ".\n"))
    }
  }

  ests <- lapply(quantiles, function(qu, which, x, margins, start, constrain=constrain, v=v, ...)
                                     mexDependence(x=x, which=which, dqu=qu, margins = margins, start=start, constrain=constrain, v=v),
                 which=which, x=x, margins = margins[[1]], start=start, constrain=constrain, v=v, referenceMargin=referenceMargin)

  boot <- lapply(ests, function(X, R, nPass, trace, ...)
                                bootmex(x=X, R=R, nPass=nPass, trace=trace, referenceMargin=referenceMargin),
                 R=R, nPass=nPass, trace=trace, referenceMargin=referenceMargin)

  res <- list(ests=ests,boot=boot,quantiles=quantiles)
  oldClass(res) <- "mexRangeFit"
  res
}

#' @export
print.mexRangeFit <- function(x, ...){
    out <- list(a = sapply(x$ests,function(x)x$dependence$coefficients[1,]),
                b = sapply(x$ests,function(x)x$dependence$coefficients[2,]))
    colnames(out$a) <- colnames(out$b) <- x$quantiles
    print(out)
    invisible(x)
}

#' @export
plot.mexRangeFit <- function(x, col=2, bootcol="grey", addNexcesses=TRUE, ...){
  ests <- x$ests
  boot <- x$boot
  quantiles <- x$quantiles
  PointEsts <- sapply(ests,function(X) coef(X$dependence))
  cof <- coef(ests[[1]]$dependence)
  whichName <- ests[[1]]$dependence$conditioningVariable
  which <- ests[[1]]$dependence$which
  data <- ests[[1]]$margins$data
  Names <- paste(rep(rownames(cof),dim(data)[2]-1),
                 paste(rep(colnames(cof),each=6),whichName,sep=" | "),sep="  ")
  R <- length(boot[[1]]$boot)

  for(i in 1:dim(PointEsts)[1]){
    if( sum((i %% 6) == 1:4) ){ # exclude plots from nuisance parameters m and s for which i mod 6 = 5,0 resp
      if(sum(PointEsts[i,])){
        Boot <- sapply(boot, function(x) sapply(x$boot, function(x) x$dependence[i]))
        ylim <- range(rbind(PointEsts[i,],Boot), na.rm=TRUE)
        plot(quantiles, PointEsts[i,], col=col, ylab=Names[i], type="b", ylim=ylim, ...)
        points(rep(quantiles,each=R),Boot,col=bootcol)
        points(quantiles, PointEsts[i,],col=col)
        if(addNexcesses){
          axis(3, at=axTicks(1), labels=sapply(axTicks(1), function(u,dat,which)sum(dat[,which] > quantile(dat[,which],u)),
                                               dat=data, which=which))
          mtext("# threshold excesses")
        }
      }
    }
  }
}


#' @export
ggplot.mexRangeFit <- function(data=NULL, mapping,
                             ylim = "auto",
                             ptcol="blue",
                             col="cornflowerblue",
                             bootcol="orange",
                             plot.=TRUE,
                             addNexcesses=TRUE,
                             textsize=4,
                             ..., environment){
    ests <- data$ests
    boot <- data$boot
    quantiles <- data$quantiles
    PointEsts <- sapply(ests,function(X) coef(X$dependence))
    cof <- coef(ests[[1]]$dependence)
    whichName <- ests[[1]]$dependence$conditioningVariable
    which <- ests[[1]]$dependence$which
    dat <- ests[[1]]$margins$data
    Names <- paste(rep(rownames(cof),dim(dat)[2]-1),
                   paste(rep(colnames(cof),each=6),whichName,sep=" | "),sep="  ")
    R <- length(boot[[1]]$boot)

    plotfn <- function(i){
        if( sum((i %% 6) == 1:4)  &  sum(PointEsts[i,])) { # exclude plots from nuisance parameters m and s for which i mod 6 = 5,0 resp
            Boot <- sapply(boot, function(x) sapply(x$boot, function(x) x$dependence[i]))

            d <- data.frame(q=quantiles,p=PointEsts[i,])
            b <- data.frame(q=rep(quantiles,each=R),b=c(Boot))

            p <- ggplot(d,aes(q,p)) +
                    labs(x="Quantiles",y=Names[i]) +
                    geom_line(colour=col) +
                    geom_point(data=b,aes(q,b),colour=bootcol,alpha=0.5) +
                    geom_point(colour=col)
        } else {
            p <- NULL
        }
        p
    }
    p <- lapply(1:dim(PointEsts)[1], plotfn)
    p <- p[!sapply(p,is.null)]

    # The loess smoother can tend to throw warnings, so suppress
    if (plot.) suppressWarnings(do.call("grid.arrange", c(p, list(ncol=2))))
    invisible(p)
}
