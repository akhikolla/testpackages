
###############################################################
#############       slopeOP and slopeSN       #################
###############################################################

#' slopeOP
#' @description Optimal partitioning algorithm for change-in-slope problem with a finite number of states (beginning and ending values of each segment is restricted to a finite set of values called states).
#' The algorithm takes into account a continuity constraint between successive segments and infers a continuous piecewise linear signal.
#' @param data vector of data to segment: a univariate time series
#' @param states vector of states = set of accessible starting/ending values for segments in increasing order.
#' @param penalty the penalty value (a non-negative real number)
#' @param constraint string defining a constraint : "null", "isotonic", "unimodal" or "smoothing"
#' @param minAngle a minimal inner angle in degree between consecutive segments in case constraint = "smoothing"
#' @param type string defining the pruning type to use. "null" = no pruning, "channel" = use monotonicity property, "pruning" = pelt-type property
#' @param testMode a boolean, if true the function also returns the percent of elements to scan (= ratio scanned elements vs. scanned elements if no pruning)
#' @return a list of 3 elements  = (changepoints, states, globalCost). (Pruning is optional)
#' \describe{
#' \item{\code{changepoints}}{is the vector of changepoints (we return the extremal values of all segments from left to right)}
#' \item{\code{states}}{is the vector of successive states. states[i] is the value we inferred at position changepoints[i]}
#' \item{\code{globalCost}}{is a number equal to the global cost of the non-penalized change-in-slope problem. That is the value of the fit to the data ignoring the penalties for adding changes}
#' \item{\emph{pruning}}{is the percent of positions to consider in cost matrix Q  (returned only if testMode = TRUE)}
#' }
#' @examples
#' myData <- slopeData(index = c(1,100,200,300), states = c(0,5,3,6), noise = 1)
#' slopeOP(data = myData, states = 0:6, penalty = 10)
slopeOP <- function(data, states, penalty = 0, constraint = "null", minAngle = 0, type = "channel", testMode = FALSE)
{
  ############
  ### STOP ###
  ############
  if(!is.numeric(data)){stop('data values are not all numeric')}
  if(!is.numeric(states)){stop('states are not all numeric')}
  if(is.unsorted(states)){stop('states must be in increasing order')}
  if(length(unique(states)) < length(states)){stop('states is not a strictly increasing sequence')}
  if(!is.double(penalty)){stop('penalty is not a double.')}
  if(penalty < 0){stop('penalty must be non-negative')}
  if(!is.double(minAngle)){stop('minAngle is not a double.')}
  if(minAngle < 0 || minAngle > 180){stop('minAngle must lie between 0 and 180')}
  allowed.constraints <- c("null", "isotonic", "unimodal", "smoothing")
  if(!constraint %in% allowed.constraints){stop('constraint must be one of: ', paste(allowed.constraints, collapse=", "))}
  allowed.types <- c("null", "channel", "pruning")
  if(!type %in% allowed.types){stop('type must be one of: ', paste(allowed.types, collapse=", "))}
  if(!is.logical(testMode)){stop('testMode must be a boolean')}

  ###CALL Rcpp functions###
  res <- slopeOPtransfer(data, states, penalty, constraint, minAngle, type)

  ###Response class slopeOP###
  ### ATTENTION : we here remove one penalty to globalCost
  if(testMode == FALSE){response <- list(changepoints = res$changepoints, parameters = res$parameters, globalCost = res$globalCost - (length(res$changepoints) - 1) * penalty)}
  if(testMode == TRUE){response <- list(changepoints = res$changepoints, parameters = res$parameters, globalCost = res$globalCost - (length(res$changepoints) - 1) * penalty, pruning = res$pruningPower)}

  attr(response, "class") <- "slopeOP"

  return(response)
}


#' slopeSN
#' @description Segment neighborhood algorithm for change-in-slope problem with a finite number of states (beginning and ending values of each segment is restricted to a finite set of values called states).
#' The algorithm takes into account a continuity constraint between successive segments and infers a continuous piecewise linear signal with a given number of segments.
#' @param data vector of data to segment: a univariate time series
#' @param states vector of states = set of accessible starting/ending values for segments in increasing order.
#' @param nbSegments the number of segments to infer
#' @param constraint string defining a constraint : "null", "isotonic"
#' @param testMode a boolean, if true the function also returns the percent of elements to scan (= ratio scanned elements vs. scanned elements if no pruning)
#' @return a list of 3 elements  = (changepoints, states, globalCost). (Pruning is optional)
#' \describe{
#' \item{\code{changepoints}}{is the vector of changepoints (we return the extremal values of all segments from left to right)}
#' \item{\code{states}}{is the vector of successive states. states[i] is the value we inferred at position changepoints[i]}
#' \item{\code{globalCost}}{is a number equal to the global cost of the non-penalized change-in-slope problem. That is the value of the fit to the data ignoring the penalties for adding changes}
#' \item{\emph{pruning}}{is the percent of positions to consider in cost matrix Q  (returned only if testMode = TRUE)}
#' }
#' @examples
#' myData <- slopeData(index = c(1,100,200,300), states = c(0,5,3,6), noise = 1)
#' slopeSN(data = myData, states = 0:6, nbSegments = 2)
slopeSN <- function(data, states, nbSegments = 1, constraint = "null", testMode = FALSE)
{
  ############
  ### STOP ###
  ############
  if(!is.numeric(data)){stop('data values are not all numeric')}
  if(!is.numeric(states)){stop('states are not all numeric')}
  if(is.unsorted(states)){stop('states must be in increasing order')}
  if(length(unique(states)) < length(states)){stop('states is not a strictly increasing sequence')}

  if(nbSegments < 1){stop('nbSegments < 1')}
  if(nbSegments%%1 > 0){stop('nbSegments is not an integer')}

  allowed.constraints <- c("null", "isotonic")
  if(!constraint %in% allowed.constraints){stop('constraint must be one of: ', paste(allowed.constraints, collapse=", "))}

  if(!is.logical(testMode)){stop('testMode must be a boolean')}
  if(nbSegments + 1 > length(data)){stop('you can not have more segments that data points')}

  ###CALL Rcpp functions###
  res <- slopeSNtransfer(data, states, nbSegments, constraint)

  ###Response class slopeOP###
  ### ATTENTION : we here remove one penalty to globalCost
  if(testMode == FALSE){response <- list(changepoints = res$changepoints, parameters = res$parameters, globalCost = res$globalCost )}
  if(testMode == TRUE){response <- list(changepoints = res$changepoints, parameters = res$parameters, globalCost = res$globalCost, pruning = res$pruningPower)}

  attr(response, "class") <- "slopeOP"

  return(response)
}


###############################################################
#############     data generator and plot     #################
###############################################################

#' slopeData
#' @description Generate data with a given continuous piecewise linear model
#' @param index a vector of increasing changepoint indices
#' @param states vector of successive states
#' @param noise noise level = standard deviation of an additional normal noise
#' @param outlierDensity probability for a datapoint to be an outlier (has to be close to 0)
#' @param outlierNoise noise level for outlier data points
#' @return a vector of simulated data
#' @examples
#' myData <- slopeData(index = c(1,100,200,300), states = c(0,5,3,6), noise = 1)
slopeData <- function(index, states, noise = 0, outlierDensity = 0, outlierNoise = 50)
{
  ############
  ### STOP ###
  ############
  if(!is.numeric(index)){stop('data values are not all numeric')}
  if(is.unsorted(index)){stop('index should be an increasing vector')}
  if(length(unique(index)) < length(index)){stop('index is not a strictly increasing sequence')}

  if(!is.numeric(states)){stop('states are not all numeric')}
  if(length(index) != length(states)){stop('index and states vectors are of different size')}
  if(!is.double(noise)){stop('noise is not a double.')}
  if(noise < 0){stop('noise must be nonnegative')}

  steps <- diff(states)/diff(index)
  response <- rep(steps, diff(index))
  response <- c(states[1], cumsum(response) + states[1])
  if(outlierDensity == 0)
  {
    response <- response + rnorm(length(response), 0, noise)
  }
  else
  {
    S <- sample(c(0,1), size = length(response), replace = TRUE, prob = c(1 - outlierDensity, outlierDensity))
    response <- response + (1-S)*rnorm(length(response), 0, noise)
    response <- response + S*rnorm(length(response), 0, outlierNoise)
  }
  return(response)
}


#' plot.slopeOP
#' @description Plot the result of the slopeOP function and the data
#' @param x a slopeOP class object
#' @param ... other parameters
#' @param data the data from which we get the slopeOP object x
#' @param chpt vector of changepoints of the model
#' @param states vector of states of the model
#' @return plot data and the inferred slopeOP result (and the model if specified in 'chpt' and 'states' parameters)
#' @examples
#' myData <- slopeData(index = c(1,100,200,300), states = c(0,5,3,6), noise = 2)
#' s <- slopeOP(data = myData, states = 0:6, penalty = 20)
#' plot(s, data = myData, chpt = c(1,100,200,300), states = c(0,5,3,6))
plot.slopeOP <- function(x, ..., data, chpt = NULL, states = NULL)
{
  n <- 1:length(data)
  p <- length(x$changepoints)
  xbis <- x$changepoints
  y <- x$parameters

  #plot the data
  plot(1:length(data), data, pch = '+')
  #plot the inferred segments IN RED
  for(i in 1:(p-1))
  {
    segments(xbis[i], y[i], xbis[i+1], y[i+1], col= 2, lty = 1, lwd = 3)
  }
  #plot the initial model IN BLUE for simulated data (see parameters in slopeData function)
  if(length(chpt) > 0 && length(chpt) == length(states))
  {
    q <- length(chpt)
    for(i in 1:(q-1))
    {
      segments(chpt[i], states[i], chpt[i+1], states[i+1], col= 4, lty = 1, lwd = 3)
    }
  }
}

