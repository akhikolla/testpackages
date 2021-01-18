#' Plot MCMC algorithm
#' 
#' @param x An Rcpp_MCMC object initialized with \code{initializeMCMCObject}.
#' 
#' @param zoom.window A vector describing the start and end of the zoom window.
#' 
#' @param what character defining if log(Posterior) (Default) or log(Likelihood) 
#' options are: LogPosterior or logLikelihood
#' 
#' @param ... Arguments to be passed to methods, such as graphical parameters.
#' 
#' @return This function has no return value.
#' 
#' @description This function will plot the logLikelihood trace, and if the Hmisc package is installed, it will 
#'  plot a subplot of the logLikelihood trace with the first few samples removed.
plot.Rcpp_MCMCAlgorithm <- function(x, what = "LogPosterior", zoom.window = NULL, ...)
{
  if(what[1] == "LogPosterior")
  {
    trace <- x$getLogPosteriorTrace()
    ylab = "log(Posterior Probability)"
  }else{
    trace <- x$getLogLikelihoodTrace()
    ylab = "log(Likelihood Probability)"
  }
  trace <- trace[-1]
  
  trace.length <- length(trace)
  
  zoomStart <- round(0.9*trace.length)
  zoomEnd <- trace.length
  logL <- mean(trace[zoomStart:trace.length])
  #TODO change main title
  plot(trace, type="l", main=paste0(ylab, ": ", logL), xlab="Sample", ylab=ylab)
  grid (NULL,NULL, lty = 6, col = "cornsilk2")
  trace[trace == -Inf] <- NA
  
  # TODO (Cedric): get rid of that line once problem with first element beeing 0 is solved
  trace <- trace[-1]
  
  if(!(is.null(zoom.window))) {
    zoomStart <- zoom.window[1]
    zoomEnd <- zoom.window[2]
  }
  else{
    warning("No window was given, zooming in at last 10% of trace")
  }
  
  Hmisc::subplot(
    plot(zoomStart:zoomEnd, trace[zoomStart:zoomEnd], type="l", xlab=NA, ylab=NA, las=2, cex.axis=0.55), 
    0.8*(round(0.9*trace.length)), (min(trace, na.rm = T)+max(trace, na.rm = T))/2, size=c(3,2))
}

#' Autocorrelation function for the likelihood or posterior trace
#' 
#' @param mcmc object of class MCMC
#' @param type "LogPosterior" or "LogLikelihood", defaults to "LogPosterior"
#' @param samples number of samples at the end of the trace used to calculate the acf
#' @param lag.max Maximum amount of lag to calculate acf. Default is 10*log10(N), where N i the number of observations.
#' @param plot logical. If TRUE (default) a plot of the acf is created
#' 
#' @description The function calculates and by defaults plots the acf and estimates the autocorrelation in the trace.
#' 
#' @seealso \code{\link{acfCSP}}
#' 
acfMCMC <- function(mcmc, type = "LogPosterior", samples = NULL, lag.max = 40, plot = TRUE)
{
  if(type == "LogPosterior")
  {
    trace <- mcmc$getLogPosteriorTrace()
  }else{
    trace <- mcmc$getLogLikelihoodTrace()
  }
  if(is.null(samples)){ samples <- round(10*log10(length(trace))) }
  
  trace <- trace[(length(trace)-samples):length(trace)]
  trace.acf <- acf(x = trace, lag.max = lag.max, plot = FALSE)
  if(plot){
    header <- paste(type, "Trace Autocorrelation",sep=" ")
    plot(x = trace.acf, xlab = "Lag time", ylab = "Autocorrelation", main = header)
  }else{
    return(trace.acf)
  }
}
