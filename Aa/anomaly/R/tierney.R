#' tierney
#' 
#' Transforms the data X by centring and scaling using \eqn{X_{ij}^{'} = \frac{X_{ij}-\mu_{ij}}{\sigma_{ij}}} where \eqn{\mu_{ij}} and \eqn{\sigma_{ij}} are robust quantile based 
#' sequential estimates for the mean and standard deviation of each variate (column) \eqn{X_{i}} of X calculated up to time j. The estimates \eqn{\mu_{ij}} and \eqn{\sigma_{ij}} are
#' calculated from sequential estimates for the median and inter-quartile range developed by Tierney et al (1983).  This method is the default value for the
#' transform argument used by the \code{\link{scapa.uv}} function.
#' 
#'
#' @param X A numeric matrix containing the data to be transformed.
#' @param burnin Specifies the period used to stabalise the quantile estimates. The default value is 10.
#' 
#' @return A numeric matrix containing the transformed data. 
#'
#' @references \insertRef{Schruben:1983:OTI:2771114.2771123}{anomaly}
#' 
#' @examples
#' library(anomaly)
#' data(machinetemp)
#' attach(machinetemp)
#' plot(temperature)
#' temperature<-tierney(temperature,burnin=4305)
#' plot(temperature)
#' @export
tierney<-function(X,burnin=10)
{
    if(is.vector(X))
    {
        ests<-sequential_ests(X,burnin)
        return((X-ests[[1]])/ests[[2]])
    }
    else if(is.matrix(X))
    {
        return(Reduce(cbind,Map(function(i) array(tierney(X[,i],burnin),c(nrow(X),1)),1:ncol(X))))  
    }
    else
    {
        # incorrect type - throw an exception
    }   
}


sequential_ests<-function(data,burnin = 10)
{
  
  # Error Traps
  ## UPDATED VERSION - LB - sequential_estimates
  # Traps for burnin
  
  if(!(is.numeric(burnin))){
    stop("burnin has to be numeric.")
  }
  
  if (length(burnin) == 0 ){
    stop("input for burnin has length 0.")
  }
  
  if (length(burnin) > 1 ){
    burnin = burnin[1]
    warning("length of input for burnin exceeds 1. Only the first entry was kept.")
  }
  
  if(!(is.numeric(burnin))){
    stop("burnin has to be numeric.")
  }
  
  if((is.infinite(burnin))){
    stop("input for burnin is infinite.")
  }
  
  if((is.nan(burnin))){
    stop("input for burnin is NaN.")
  }
  
  if(!(is.integer(burnin))){
    burnin = as.integer(burnin)
    # warning("non-integer input for burnin. The input was converted to an integer using as.integer.")
  }
  
  if (burnin < 10){
    stop("argument for burnin needs to be at least 10.")
  }
  
  # Traps for data
  
  if(!(is.numeric(data))){
    stop("data has to be numeric.")
  }
  
  data = as.vector(data)
  
  if(length(data) <= burnin){
    stop("length of data less than burnin.")
  }
  
  if(sum(is.na(data)) > 0){
    stop("Input for data contains NAs.")
  }
  
  if(sum(is.nan(data)) > 0){
    stop("Input for data contains NaNs.")
  }
  
  if(sum(is.infinite(data)) > 0){
    stop("Input for data contains Infinite entries.")
  }
  
  
  # Return 
  
  # initial estimates done in R
  slq = quantile(data[1:burnin], probs = 0.25)
  smed = quantile(data[1:burnin], probs = 0.5)
  suq = quantile(data[1:burnin], probs = 0.75)
  
  scale = IQR(data[1:burnin])
  c = (scale/burnin) * sum( (1:burnin)^(-0.5) )
  
  flq = (1/(2*c*burnin)) * max( sum( abs(data[1:burnin] - slq) <= c ), 1)
  fmed = (1/(2*c*burnin)) * max( sum( abs(data[1:burnin] - smed) <= c ), 1)
  fuq = (1/(2*c*burnin)) * max( sum( abs(data[1:burnin] - suq) <= c ), 1)
  
  n = length(data)
  return(marshall_sequential_ests(data, n, burnin, slq, flq, smed, fmed, suq, fuq))
  
} 
