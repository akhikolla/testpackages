#' Estimate the EGP3 distribution power parameter over a range of thresholds
#'
#' Estimate extended generalized Pareto distribution power parameter over a
#' range of values, using maximum (penalized) likelihood.
#'
#' Papastathopoulos and Tawn present 3 extended versions of the generalized
#' Pareto distribution. Using the \code{egp3} texmex family object, the power
#' parameter in the EGP3 distribution is estimated on the log scale, a
#' confidence interval is calculated and the result is transformed back to the
#' scale of the power parameter and returned to the user.
#'
#' When the power paramer, kappa, is equal to 1, the EPG3 distribution is
#' identical to the generalized Pareto distribution. Therefore, the plot of the
#' estimated parameter over a range of thresholds provides a diagnostic for
#' threshold selection: the lowest value of kappa whose confidence interval
#' includes 1 is suggested as the threshold for generalized Pareto modelling.
#'
#' If lower thresholds are used and the EGP3 distribution itself is used for
#' modelling, some care should be taken to ensure the model provides a
#' reasonable degree of fit to the data. Limited experience suggests that such
#' models seldom fit well and the main value of the EGP3 distribution is as a
#' diagnostic for threshold selection as described here.
#'
#' Note this function does not extend to assessing model fit when there are
#' covariates included in the model.
#'
#' @aliases egp3RangeFit print.egp3RangeFit plot.egp3RangeFit ggplot.egp3RangeFit
#' @usage egp3RangeFit(data, umin=quantile(data, .05), umax=quantile(data,
#' .95), nint = 10, penalty = "gaussian", priorParameters = NULL, alpha=0.05)
#' \method{print}{egp3RangeFit}(x, ...)
#' \method{plot}{egp3RangeFit}(x, xlab = "Threshold", ylab = "kappa", main = NULL, addNexcesses=TRUE, log.="", ...)
#' \method{ggplot}{egp3RangeFit}(data, mapping, xlab = "Threshold", ylab = expression(kappa),
#' main=NULL,fill="orange", col="blue",addNexcesses=TRUE, textsize=4, ..., environment)
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
#' @param x Argument to the \code{print} functions.
#' @param xlab Label for the x-axis.
#' @param ylab Label for the y-axis.
#' @param main The main title.
#' @param textsize Size of text for annotation showing number of threshold excesses.
#' @param addNexcesses Annotate top axis with numbers of threshold excesses
#' arising with the corresponding values of threshold on the bottom axis.
#' @param log. Argument passed through to \code{plot}. Can take values "x" for
#' plotting the x-axis on the log scale, "y" for plotting the y-axis on the log
#' scale, "xy" for both, or "" (the default) for neither.
#' @param mapping,fill,col,environment Arguments to ggplot method.
#' @param \dots Arguments to \code{plot}.
#' @author Harry Southworth
#' @seealso \code{\link{evm}}, \code{\link{gpdRangeFit}}, \code{\link{mrl}}
#' @references I. Papastathopoulos and J. A. Tawn, Extended generalized Pareto
#' modles for tail estimation, Journal of Statistical Planning and Inference,
#' 143, 131 -- 143, 2013
#' @keywords models
#' @examples
#'
#' \donttest{ # because of the time it takes to run
#' erf <- egp3RangeFit(rain)
#' plot(erf)
#' ggplot(erf)
#' }
#'
#' @export egp3RangeFit
egp3RangeFit <-
function (data, umin=quantile(data, .05), umax=quantile(data, .95),
          nint = 10,
          penalty="gaussian", priorParameters=NULL, alpha=.05) {

  #if (umin < 0)
  #  stop("umin < 0: data must be non-negative. Add a constant or increase umin and try again")

  m <- s <- hi <- lo <- rep(0, nint)
  u <- seq(umin, umax, length = nint)
  qz <- qnorm(1 - alpha/2)

  for (i in 1:nint) {
    z <- evm(data, th=u[i], penalty=penalty, priorParameters=priorParameters, family=egp3)
    m[i] <- z$coefficients[1]
    s[i] <- z$se[1]
  }

  # egp3 family works with labmda = log(kappa)
  hi <- exp(m + qz * s)
  lo <- exp(m - qz * s)

  res <- list(th=u, par=exp(m) , hi=hi, lo=lo, data=data)
  oldClass(res) <- 'egp3RangeFit'
  res
}

#' @export
print.egp3RangeFit <- function(x, ...){
  print(cbind(threshold=x$th, kappa=x$par, lo=x$lo, hi=x$hi))
  invisible(x)
}

#' @export
plot.egp3RangeFit <- function(x, xlab="Threshold", ylab="kappa",
                              main=NULL, addNexcesses=TRUE, log.="", ...){

  yl <- range(x$hi, x$lo, 1)
  plot(x$th, x$par, ylim = yl, type = "b",
       xlab=xlab, ylab=ylab, main=main, log=log., ...)
  for (j in 1:length(x$th)){
    lines(c(x$th[j], x$th[j]), c(x$hi[j], x$lo[j]))
  }
  if(addNexcesses){
    axis(3, at=axTicks(1), labels=sapply(axTicks(1), function(u) sum(x$data > u)), cex=0.5)
    mtext("# threshold excesses")
  }
  abline(h=1, lty=2)

  invisible()
}

#' @export
ggplot.egp3RangeFit <- function(data, mapping, xlab = "Threshold", ylab = expression(kappa), main=NULL,
                                fill="orange", col="blue",
                                addNexcesses=TRUE, textsize=4, ..., environment)
{
    d <- data.frame(hi=data$hi,lo=data$lo,th=data$th,par=data$par)
    poly <- data.frame(x=c(data$th, rev(data$th)), y=c(data$lo, rev(data$hi)))

    p <- ggplot(data=d,aes(th,par)) +
        geom_polygon(data=poly,aes(x,y),fill=fill, alpha=.5) +
        geom_line(colour=col) +
        geom_hline(yintercept=1,linetype=2) +
        labs(x=xlab,y=ylab,title=main)

    if (addNexcesses)
        p <- addExcesses(p, poly$x, poly$y, data=data$data, textsize=textsize)

    invisible(p)
}
