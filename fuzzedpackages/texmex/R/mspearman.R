#---------------------------------------------------------------------------
#Function MCS() evaluates the MCS wrt to Schmid et al. 2006 equation (17)
#First and second args are
#X (matrix) and p (vector of probabilities).
#Matrix X has to be i*j where i=1,...,n, and j=1,...,d and d,n
#denote dimension and length, respectively.	
#---------------------------------------------------------------------------

.MCSlower <- function(U, p)
  {
    d <- dim(U)[1]
    n <- dim(U)[2]
    res1         <- p - U
    res1[res1 < 0] <- 0
    res2         <- apply(res1,2,prod)
    res3         <- sum(res2)
    res4         <- ( (1/n)*res3-((p^2)/2)^(d) )/
      ( (p^(d+1))/(d+1) -((p^2)/2)^(d))
    return(res4)
  }


#' Multivariate conditional Spearman's rho
#' 
#' Compute multivariate conditional Spearman's rho over a range of quantiles.
#' 
#' The method is described in detail by Schmid and Schmidt (2007).  The main
#' code was written by Yiannis Papastathopoulos, wrappers written by Harry
#' Southworth.
#' 
#' When the result of a call to \code{bootMCS} is plotted, simple quantile
#' bootstrap confidence intervals are displayed.
#' 
#' @aliases MCS plot.MCS ggplot.MCS print.MCS bootMCS plot.bootMCS print.bootMCS summary.bootMCS 
#' 
#' @param X A matrix of numeric variables.
#' @param p The quantiles at which to evaluate.
#' @param R The number of bootstrap samples to run. Defaults to \code{R = 100}.
#' @param trace How often to inform the user of progress. Defaults to
#' \code{trace = 10}.
#' @param x,object An object of class \code{MCS} or \code{bootMCS}.
#' @param xlab,ylab Axis labels.
#' @param alpha A 100(1 - alpha)\% pointwise confidence interval will be
#' produced.  Defaults to \code{alpha = 0.05}.
#' @param ylim Plotting limits for bootstrap plot.
#' @param data,mapping,main,environment Arguments to ggplot method.
#' @param ... Optional arguments to be passed into methods.
#' @return MCS returns an object of class \code{MCS}.  There are plot and
#' print methods available for this class.
#' 
#' \item{MCS }{The estimated correlations.} \item{p }{The quantiles at which
#' the correlations were evaluated at} \item{call}{The function call used.}
#' 
#' bootMCS returns an object of class \code{bootMCS}. There are plot and
#' summary methods available for this class.
#' 
#' \item{replicates}{Bootstrap replicates.} \item{p }{The quantiles at which
#' the correlations were evaluated at} \item{R}{Number of bootstrap samples.}
#' \item{call}{The function call used.}
#' 
#' @author Yiannis Papastathopoulos, Harry Southworth
#' @seealso \code{\link{chi}}
#' @references F. Schmid and R. Schmidt, Multivariate conditional versions of
#' Spearman's rho and related measures of tail dependence, Journal of
#' Multivariate Analysis, 98, 1123 -- 1140, 2007
#' @keywords multivariate
#' @examples
#' 
#' D <- liver[liver$dose == "D",]
#' plot(D)
#' \donttest{
#' Dmcs <- bootMCS(D[, 5:6])
#' Dmcs
#' plot(Dmcs)
#' }
#' @export MCS
MCS <- function(X,p=seq(.1, .9, by=.1)) {
    theCall <- match.call()
     
    U <- t(apply(X, 2, edf)) #transpose cause apply transposes g(X), g:edf
    n    <- length(p)
    res <- vapply(p, FUN=.MCSlower, FUN.VALUE=0, U=U)

    res <- list(mcs=res, p=p, call=theCall)
    oldClass(res) <- "MCS"
    res
  }

#' @rdname MCS
#' @export
plot.MCS <- function(x, xlab="p", ylab= "MCS", ...){
   plot(x$p, x$mcs, type="l", xlab=xlab, ylab=ylab, ...)
   invisible()
}

#' @rdname MCS
#' @export
ggplot.MCS <- function(data, mapping, main="", ..., environment){
   ggplot(data=data.frame(p=data$p, MCS=data$mcs), aes(p,MCS)) + geom_line() + ggtitle(main)
}

#' @export
print.MCS <- function(x, ...){
    print(x$call)
    cat("Multivariate conditional Spearman's rho.\n\n", sep = "")
    res <- x$mcs
    names(res) <- x$p
    print(res)
    invisible(x)
}

#------------------------------------------------
#Bootstrap
#------------------------------------------------
#' @rdname MCS
#' @export
bootMCS <- function(X,p=seq(.1, .9, by=.1),R=100, trace=10) {
   theCall <- match.call()
   bfun <- function(i, data, p, trace){
       if (i %% trace == 0){ message("Replicate ", i) }
       d <- data[sample(1:nrow(data), replace=TRUE),]
       MCS(d, p)$mcs
   }

   res <- sapply(1:R, bfun, data=X, p=p, trace=trace)
   res <- list(replicates=res, p=p, R=R, call=theCall)
   oldClass(res) <- "bootMCS"
   invisible(res)
}

#' @rdname MCS
#' @export
ggplot.bootMCS <- function(data, mapping, main="", alpha=.05, ylim, ..., environment){
    ci <- apply(data$replicates, 1, quantile, prob=c(1-alpha/2, alpha/2))

    poly <- data.frame(p=c(data$p,rev(data$p)),
                       ci = c(ci[1,],rev(ci[2,])))
    dat <- data.frame(p=data$p,MCS = rowMeans(data$replicates))

    if (missing(ylim)) ylim <- range(ci)
    
    p <- ggplot(dat,aes(p,MCS)) + geom_line(colour="blue") + 
        geom_polygon(data=poly,mapping=aes(p,ci),fill="orange",alpha=0.5) + 
        coord_cartesian(ylim=ylim) +
        labs(title=main,
             x=paste("p\n",100*(1-alpha), "% interval. ", data$R, " bootstrap samples were performed", sep=""))

    p
}

#' @rdname MCS
#' @export
plot.bootMCS <- function(x, xlab="p", ylab= "MCS",alpha=.05, ylim, ...){
   m <- rowMeans(x$replicates)
   ci <- apply(x$replicates, 1, quantile, prob=c(1-alpha/2, alpha/2))
   if (missing(ylim)){ ylim <- range(ci) }
   plot(x$p, m, type="l", ylim=ylim,
        xlab=xlab, ylab=ylab,
        sub=paste(100*(1-alpha), "% interval. ", x$R, " bootstrap samples were performed", sep=""),...)
   lines(x$p, ci[1,], lty=2)
   lines(x$p, ci[2,], lty=2)
   invisible(ci)
}

#' @export
print.bootMCS <- function(x, ...){
    print(x$call)
    cat("Multivariate conditional Spearman's rho.\n", x$R, " bootstrap samples were performed.\n\n",
        sep = "")
     
    m <- rowMeans(x$replicates)
    names(m) <- x$p
    print(m)
    invisible(x)
}

#' @rdname MCS
#' @export
summary.bootMCS <- function(object, alpha=.05, ...){
    m <- rowMeans(object$replicates)
    ci <- apply(object$replicates, 1, quantile, prob=c(alpha/2, 1 - alpha/2))
    res <- cbind(m, t(ci))
    dimnames(res) <- list(object$p, c("Mean", rownames(ci)))
    res <- list(res=res,R=object$R)
    oldClass(res) <- "summary.bootMCS"
    res
}

#' @rdname MCS
#' @export
print.summary.bootMCS <- function(x, ...){
    cat("Multivariate conditional Spearman's rho.\n", x$R, " bootstrap samples were performed.\n\n",
        sep = "")
    print(x$res)
    invisible(x)
}

