#' Extremal index estimation and automatic declustering
#'
#' Given a threshold which defines excesses above that threshold, estimate the
#' extremal index of a dependent sequence by using the method of Ferro and
#' Segers, 2003.  The extremal index estimate can then be used to carry out
#' automatic declustering of the sequence to identify independent clusters and
#' estimate the GPD for cluster maxima.  Graphical diagnostics of model fit are
#' available.
#'
#' The function \code{extremalIndex} estimates the extremal index of a
#' dependent series of observations above a given threshold \code{threshold},
#' returning an object of class "extremalIndex".  Plot and print methods are
#' available for this class. A graphical diagnostic akin to Figure 1 in Ferro
#' and Segers (2003) is produced by the \code{plot} method for this class.
#' This plot is used to test the model assumption underpinning the estimation,
#' with good fit being indicated by interexceedance times which correspond to
#' inter-cluster times lying close to the diagonal line indicated.
#'
#' In addition to good model fit, an appropriate choice of threshold is one
#' above which the estimated extremal index is stable over further, higher
#' thresholds (up to estimation uncertainty).  This can be assessed by using
#' the function \code{extremalIndexRangeFit}, which examines a range of
#' threshold values.  At each threshold, the extremal index is estimated; that
#' estimate is used to decluster the series and the parameters of the GPD are
#' optionally estimated for the resulting declustered series.  Uncertainty in
#' the estimation of the extremal index and GPD parameters is assessed by using
#' a bootstrap scheme which accounts for uncertainty in the extremal index
#' estimation, and the corresponding uncertainty in the declustering of the
#' series. There are \code{plot} and \code{ggplot} methods for output of this function, which is of class \code{extremalIndexRangeFit}.
#'
#' The function \code{declust} returns an object of class "declustered",
#' identifying independent clusters in the original series. Print, plot and
#' show methods are available for this class. The GPD model can be fitted to
#' objects of this class, including the use of covariates in the linear
#' predictors for the parameters of the GPD.  See examples below.
#'
#' @aliases extremalIndex extremalIndexRangeFit declust declust.default
#' declust.extremalIndex print.extremalIndex plot.declustered print.declustered bootExtremalIndex evm.declustered plot.extremalIndexRangeFit ggplot.extremalIndexRangeFit
#' @usage extremalIndex(y, data = NULL, threshold)
#'
#' extremalIndexRangeFit(y, data = NULL, umin = quantile(y,.5), umax =
#' quantile(y, 0.95), nint = 10, nboot = 100, alpha = .05, estGPD=TRUE,
#' verbose = TRUE, trace = 10, ...)
#'
#' bootExtremalIndex(x)
#'
#' declust(y, r=NULL, data = NULL, ...)
#'
#' \method{declust}{extremalIndex}(y, r=NULL,...)
#'
#' \method{plot}{declustered}(x, ylab = "Data",...)
#'
#' \method{evm}{declustered}(y, data=NULL, family=gpd, ...)
#'
#' \method{plot}{extremalIndexRangeFit}(x,addNexcesses=TRUE,estGPD=TRUE,...)
#'
#' \method{print}{extremalIndex}(x,...)
#'
#' \method{print}{declustered}(x,...)
#'
#' \method{ggplot}{extremalIndexRangeFit}(data=NULL, mapping, xlab, ylab, main,
#' ylim = "auto",ptcol="dark blue",col="dark blue",fill="orange",
#' textsize=4,addNexcesses=TRUE,estGPD=TRUE,..., environment)
#'
#' @param y Argument to function \code{extremalIndex}: either a numeric vector
#' or the name of a variable in \code{data}.
#' @param data A data frame containing \code{y} and any covariates. In
#' \code{evm.declustered}, it should be NULL and is included to match the
#' arguments of generic \code{evm}.
#' @param family The type of extreme value model. The user should not change
#' this from its default in \code{evm.declustered}.
#' @param threshold The threshold for \code{y}, exceedances above which will be
#' used to estimate the extremal index and carry out automatic declustering.
#' @param x Objects passed to methods.
#' @param r Positivie integer: run length to be used under "runs" declustering.
#' If specified then so-called "runs" declustering will be carried out,
#' otherwise defaults to NULL in which case the automatic "intervals"
#' declustering method of Ferro and Segers is used.
#' @param umin The minimum threshold above which to esimate the parameters.
#' @param umax The maximum threshold above which to esimate the parameters.
#' @param nint The number of thresholds at which to perform the estimation.
#' @param nboot Number of bootstrap samples to simulate at each threshold for
#' estimation.
#' @param alpha 100(1 - alpha)\% confidence intervals will be plotted with the
#' point estimates. Defaults to \code{alpha = 0.05}.
#' @param xlab Label for the x-axis (ggplot).
#' @param ylab Label for the y-axis (ggplot).
#' @param addNexcesses Whether to annotate the top axis of plots with the
#' number of excesses above the corresponding threhsold. Defaults to
#' \code{TRUE}.
#' @param estGPD Whether to estimate GPD parameters at each choice of
#' thereshold -- defaults to \code{TRUE} in which case the GPD parameters are
#' estimated.
#' @param verbose Whether to report on progress in RangeFit calculations.
#' Defaults to \code{TRUE}.
#' @param trace How frequently to report bootstrap progress in RangeFit
#' calculations.  Defaults to 10.
#' @param  mapping,main,ylim,ptcol,col,fill,textsize,environment Further arguments to ggplot method.
#' @param ... Further arguments to methods.
#' @return The function \code{extremalIndex} returns a list of class
#' "extremalIndex": \item{EIintervals}{Estimate of the extremal index by using
#' the intervals estimator of Ferro and Segers.} \item{threshold}{threshold for
#' declustering and estimation} \item{TotalN}{length of original data series}
#' \item{nExceed}{number of exceedances of \code{threshold} in original
#' series.} \item{thExceedanceProb}{probablity of threshold exceedance in
#' original series.} \item{call}{the original function call }
#' \item{interExceedTimes}{times between threshold exceedances}
#' \item{thExceedances}{observation from the original series which are above
#' \code{threshold}} \item{exceedanceTimes}{times of occurrance of threshold
#' exceedances} \item{y}{original dependent series} \item{data}{data frame or
#' NULL}
#'
#' The function \code{declust} returns a list of type "declustered":
#'
#' \item{clusters}{integer labels assigning threshold exceedances to clusters}
#' \item{sizes}{number of exceedances in each cluster}
#' \item{clusterMaxima}{vector made up of the largest observation from each
#' distinct cluster.  In the case of ties, the first value is taken.}
#' \item{isClusterMax}{logical; length equal to number of threshold
#' exceedances, value is \code{TRUE} for threshold exceedances which correspond
#' to cluster maxima} \item{y}{see entry for object of class "extremalIndex"
#' above} \item{data}{see entry for object of class "extremalIndex" above}
#' \item{threshold}{see entry for object of class "extremalIndex" above}
#' \item{EIintervals}{see entry for object of class "extremalIndex" above}
#' \item{call}{see entry for object of class "extremalIndex" above}
#' \item{InterExceedTimes}{times between threshold exceedances, length is one
#' less than the number of threshold exceedances} \item{InterCluster}{logical:
#' indicates inter exceedance times larger than \code{r} the run length used
#' for declustering} \item{thExceedances}{see entry for object of class
#' "extremalIndex" above} \item{exceedanceTimes}{see entry for object of class
#' "extremalIndex" above} \item{r}{run length used for declustering}
#' \item{nClusters}{Number of indenendent clusters identified}
#' \item{method}{Method used for declustering (either "intervals" or "runs")}
#'
#' The function \code{bootExtremalIndex} return a single vector corersponding
#' to a bootstrap sample from the original series: observations are censored at
#' \code{threshold} so that values below this threshold are indicated by the
#' value -1.
#'
#' The method \code{evm} for class "declustered" returns an object of type
#' "evmOpt" or "evmSim" depending on the precise function call - see
#' documentation for \code{\link{evm}}.
#' @author Janet E. Heffernan
#' @seealso \code{\link{evm}}
#' @references Ferro, C.A.T. and Segers, J., (2003) "Inference for clusters of
#' Extreme Values", JRSS B 65, Part 2, pp 545--556.
#' @examples
#'
#' par(mfrow=c(2,2));
#' extremalIndexRangeFit(summer$O3,nboot=10)
#' ei <- extremalIndex(summer$O3,threshold=45)
#' plot(ei)
#' d <- declust(ei)
#' plot(d)
#' evm(d)
#'
#' ## fitting with covariates:
#'
#' so2 <- extremalIndex(SO2,data=winter,threshold=15)
#' plot(so2)
#' so2 <- extremalIndex(SO2,data=winter,threshold=20)
#' plot(so2) ## fits better
#'
#' so2.d <- declust(so2)
#' par(mfrow=c(1,1)); plot(so2.d)
#' so2.d.gpd <- evm(so2.d) # AIC 661.1
#'
#' evm(so2.d,phi=~NO)
#' evm(so2.d,phi=~NO2)
#' evm(so2.d,phi=~O3) # better AIC 651.9
#' evm(so2.d,phi=~PM10)
#'
#' so2.d.gpd.o3 <- evm(so2.d,phi=~O3)
#'
#' par(mfrow=c(2,2)); plot(so2.d.gpd.o3)
#'
#' @export extremalIndex
extremalIndex <- function(y,data=NULL,threshold)
# intevals estimator of the Extremal Index, Ferro and Segers JRSS B (2003)
# assumes data points equally spaced in time and no missing data (ie missing time points)
{
  if (!missing(data)) {
     y <- ifelse(deparse(substitute(y))== "substitute(y)", deparse(y),deparse(substitute(y)))
     y <- formula(paste(y, "~ 1"))
     y <- model.response(model.frame(y, data=data))
  }
  theCall <- match.call()
  timeAll <- 1:length(y)
  thExceedance <- y > threshold
  thExceedanceProb <- mean(thExceedance)
  nExceed <- sum(thExceedance)
  exceedanceTimes <- timeAll[thExceedance]
  interExceedTimes <- exceedanceTimes[-1] - exceedanceTimes[-nExceed]

  if ( max(interExceedTimes) > 2 ) {
    Est <- (2 * sum(interExceedTimes - 1)^2) / ( (nExceed-1) * sum( (interExceedTimes - 1) * (interExceedTimes-2) ) )
  } else {
    Est <- (2 * sum(interExceedTimes)^2) / ( (nExceed-1) * sum( interExceedTimes^2) )
  }

  EIintervals <- min(1,Est)

  res <- list(EIintervals = EIintervals,
              threshold = threshold,
              TotalN = length(y),
              nExceed = nExceed,
              thExceedanceProb = thExceedanceProb,
              call = theCall,
              interExceedTimes = interExceedTimes,
              thExceedances = y[thExceedance],
              exceedanceTimes = timeAll[thExceedance],
              y = y,
              data = data)

  oldClass(res) <- "extremalIndex"

  res
}

#' @export
print.extremalIndex <- function(x,...)
{
  cat("\nLength of original series",x$TotalN,"\n")
  cat("Threshold", x$threshold,"\n")
  cat("Number of Threshold Exceedances",x$nExceed,"\n")
  cat("Intervals estimator of Extremal Index", x$EIintervals,"\n")
  invisible(x)
}

#' @export
plot.extremalIndex <- function(x,...)
{
  NormInterExceedTimes <- x$interExceedTimes * x$thExceedanceProb

  StdExpQuantiles <- qexp(ppoints(NormInterExceedTimes))
  Theta <- x$EIintervals

  plot(StdExpQuantiles, sort(NormInterExceedTimes),xlab="Standard Exponential Quantiles",ylab="Interexceedance Times",cex=0.7,...)
  abline(v=qexp(1-Theta))
  abline(a = -qexp(1-Theta)/Theta, b=1/Theta)
  title(paste("Threshold=",x$threshold))
  invisible()
}

#' @export
declust <- function(y, r=NULL, data=NULL, ...)
{
  if (!missing(data)) {
     y <- deparse(substitute(y))
     y <- formula(paste(y, "~ 1"))
     y <- model.response(model.frame(y, data=data))
  }
  UseMethod("declust", y)
}

#' @export
declust.default <- function(y,r=NULL,data=NULL,verbose=TRUE,...)
{
  if(missing(data)){
    ei <- extremalIndex(y,...)
  } else {
    ei <- extremalIndex(substitute(y),data,...)
  }

  declust(ei, r=r)
}

#' @export
declust.extremalIndex <- function(y,r=NULL,...)
{
  theCall <- match.call()
  Times <- y$interExceedTimes
  sTimes <- sort(Times, decreasing=TRUE)

  if(is.null(r)){
    C <- floor(y$EIintervals * y$nExceed) + 1
    C <- min(C,length(Times)) # needed if short series and C < number of interexceedance times
    while(sTimes[C-1] == sTimes[C]) C <- C-1
    r <- sTimes[C]
    method <- "intervals"
  } else {
    method <- "runs"
  }

  clusters <- rep(1,length(y$thExceedances))
  clusters[-1] <- 1+cumsum(Times > r)
  sizes <- tabulate(clusters)
  C <- max(clusters)

  clusterMaxima <- sapply(1:C,function(i) max(y$thExceedances[clusters == i]))
  isClusterMax <- rep(FALSE,length(clusters))
  for(i in 1:C){
    isClusterMax[clusters == i & y$thExceedances == max(y$thExceedances[clusters == i])][1] <- TRUE
  }

  res <- list(clusters = clusters,
              sizes=sizes,
              clusterMaxima = clusterMaxima,
              isClusterMax = isClusterMax,
              y = y$y,
              data = y$data,
              threshold=y$threshold,
              EIintervals = y$EIintervals,
              call=theCall,
              InterExceedTimes=Times,
              InterCluster = Times > sTimes[C],
              thExceedances = y$thExceedances,
              exceedanceTimes = y$exceedanceTimes,
              r=r, nClusters = C, method=method)

  oldClass(res) <- "declustered"

  res
}

#' @export
print.declustered <- function(x,...){
  print(x$call)
  cat("\nThreshold ",x$threshold,"\n")
  cat("Declustering using the",x$method,"method, run length",x$r,"\n")
  cat("Identified",length(x$sizes),"clusters.\n")
  invisible(x)
}

#' @export
plot.declustered <- function(x,ylab="Data",...){
  plot(x$y,xlab="",ylab=ylab)
  abline(h=x$threshold,col=2)
  for(i in 1:length(x$sizes)){
    points(x$exceedanceTimes[x$clusters == i],x$thExceedances[x$clusters == i],col=2,type="b")
  }
}

#' @export
bootExtremalIndex <- function(x){
  if( class(x) == "extremalIndex"){
    x <- declust(x)
  } else if(class(x) != "declustered"){
    stop("x must be of class extremalIndex or declust")
  }
  nc <- length(x$sizes)
  boot.interExceedTimes <- boot.thExceedances <- NULL
  for(i in 1:nc){
    clust <- sample(unique(x$clusters),1)
    boot.interExceedTimes <- c(boot.interExceedTimes,diff(x$exceedanceTimes[x$clusters == clust]))
    boot.thExceedances <- c(boot.thExceedances,x$thExceedances[x$clusters == clust])
    if(i < nc){
      boot.interExceedTimes <- c(boot.interExceedTimes,sample(x$InterExceedTimes[x$InterCluster], 1))
    }
  }

  boot.exceedanceTimes <- cumsum(c(1,boot.interExceedTimes))
  boot.data <- rep(-1,max(boot.exceedanceTimes))
  boot.data[boot.exceedanceTimes] <- boot.thExceedances

  boot.data
}

#' @export
extremalIndexRangeFit <- function(y,data=NULL,umin=quantile(y,.5),umax=quantile(y,0.95),nint=10,nboot=100,alpha=.05, estGPD=TRUE, verbose=TRUE, trace=10, ...){

  if (!missing(data)) {
     y <- deparse(substitute(y))
     y <- formula(paste(y, "~ 1"))
     y <- model.response(model.frame(y, data=data))
  }

  EI <- SH <- SC <-
    list(m=numeric(nint),boot=matrix(0,nrow=nint,ncol=nboot))

  u <- seq(umin, umax, length = nint)
  for (i in 1:nint) {
    if(verbose){
      cat("\n", i,"of",nint,"thresholds: bootstrap ... ")
    }
    z <- extremalIndex(y,threshold=u[i])
    EI$m[i] <- z$EIintervals
    d <- declust(z)
    if(estGPD){
      gpd.d <- evm.declustered(d)
      co.d <- coef(gpd.d)
      SH$m[i] <- co.d[2]
      SC$m[i] <- exp(co.d[1]) - co.d[2]*u[i]
    }

    for(j in 1:nboot){
      if(verbose & j %% trace == 0){
        cat(j,"")
      }
      boot <- bootExtremalIndex(d)
      z.b <- extremalIndex(boot,threshold=u[i])
      EI$boot[i,j] <- z.b$EIintervals
      if(estGPD){
        z.d <- declust(z.b)
        z.d$clusterMaxima <- rgpd(z.d$nClusters,exp(co.d[1]),co.d[2],u=z.d$threshold)
        gpd.b <- try(evm.declustered(z.d,cov="numeric"))
        if(class(gpd.b) == "try-error"){
          SH$boot[i,j] <- SC$boot[i,j] <- NA
        } else {
          SH$boot[i,j] <- coef(gpd.b)[2]
          SC$boot[i,j] <- exp(coef(gpd.b)[1]) -  SH$boot[i,j]*u[i]
        }
      }
    }
  }
  EI$ul <- apply(EI$boot,1,quantile,alpha/2)
  EI$up <- apply(EI$boot,1,quantile,1-alpha/2)
  if(estGPD){
    SC$ul <- apply(SC$boot,1,quantile,alpha/2,na.rm=TRUE)
    SC$up <- apply(SC$boot,1,quantile,1-alpha/2,na.rm=TRUE)
    SH$ul <- apply(SH$boot,1,quantile,alpha/2,na.rm=TRUE)
    SH$up <- apply(SH$boot,1,quantile,1-alpha/2,na.rm=TRUE)
  }

  EI$u <- SC$u <- SH$u <- u
  res <- list(EI=EI,SC=SC,SH=SH,y=y)
  oldClass(res) <- "extremalIndexRangeFit"
  invisible(res)
}

#' @export
ggplot.extremalIndexRangeFit <- function(data=NULL, mapping, xlab, ylab, main,
                                         ylim = "auto",
                                         ptcol="dark blue",
                                         col="dark blue",
                                         fill="orange", textsize=4,
                                         addNexcesses=TRUE,estGPD=TRUE,
                                         ..., environment){
    plots <- function(l,y,main,xlab,ylab,...){
        data <- data.frame(u=l$u,m=l$m,ul=l$ul,u=l$up)
        p <- ggplot(data,aes(u,m)) + geom_point(colour=ptcol) + labs(x=xlab,y=ylab,title=main)

        for (i in 1:dim(data)[1]){
            d <- data.frame(x=rep(l$u[i],2),y=c(l$ul[i],l$up[i]))
            p <- p + geom_line(data=d,aes(x,y),colour=col)
        }
        if (addNexcesses){
            p <- addExcesses(p, l$u, c(l$ul,l$up), data=y, textsize=textsize)
        }
    }


    res <- list(p1 = plots(data$EI,data$y,main="Extremal Index",xlab="Threshold",ylab=expression(theta),...))
    if(estGPD){
        res$p2 <- plots(data$SC,data$y,main="Scale parameter",xlab="Threshold",ylab=expression(sigma),...)
        res$p3 <- plots(data$SH,data$y,main="Shape parameter",xlab="Threshold",ylab=expression(xi),...)
    }

    invisible(res)
}

#' @export
plot.extremalIndexRangeFit <- function(x,addNexcesses=TRUE,estGPD=TRUE,...){
    plots <- function(l,y,...){
        plot(l$u, l$m, ylim=c(min(l$ul),max(l$up)),type = "b", ...)
        for (i in 1:length(l$u)) lines(c(l$u[i], l$u[i]), c(l$ul[i], l$up[i]))
        if(addNexcesses){
            axis(3,at=axTicks(1),labels=sapply(axTicks(1),function(u)max(declust(extremalIndex(y,threshold=u))$clusters)))
            mtext("# threshold excesses")
        }
    }

    plots(l=x$EI,y=x$y,main="Extremal Index",xlab="Threshold",ylab=expression(theta),...)
    if(estGPD){
        plots(l=x$SC,y=x$y,main="Scale parameter",xlab="Threshold",ylab=expression(sigma),...)
        plots(l=x$SH,y=x$y,main="Shape parameter",xlab="Threshold",ylab=expression(xi),...)
    }

}

#' @export
evm.declustered <- function(y, data=NULL, family=gpd, ...){
  myCall <- match.call()

  if(is.null(y$data)){
    res <- evm(y$clusterMaxima, th = y$threshold, family=family, ...)
  } else {
    response <- y$clusterMaxima
    dat <- cbind(response,y$data[y$y>y$threshold,][y$isClusterMax,])
    res <- evm(response, data=dat, th = y$threshold, family=family, ...)
  }

  clusterRate <- max(y$clusters) / length(y$y)
  if(class(res) == "evmOpt"){
    res$rate <- clusterRate
  } else if(class(res) == "evmSim") {
    res$map$rate <- clusterRate
  }
  res$call <- myCall
  res
}

