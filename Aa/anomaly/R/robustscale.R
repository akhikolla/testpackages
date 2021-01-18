#' robustscale
#'
#' Transforms the data X by centring and scaling using \eqn{X_{ij}^{'} = \frac{X_{i}-\mu_{i}}{\sigma_{i}}} where \eqn{\mu_{i}} and \eqn{\sigma_{i}} are robust estimates for
#' the mean and standard deviation of each variate (column), \eqn{X_{i}}, of the multivariate time series X. The estimates are calculated using the median and median absolute deviation. 
#' This method is the default value for the
#' transform argument used by the \code{\link{capa}} function, since the capa method assumes that the typical distribution of the data is standard normal.
#' 
#'
#' @param X A numeric matrix containing the data to be transformed. Each column corresponds to a component and each row to an observation.
#' 
#' @return A numeric matrix containing the transformed data. 
#' 
#' @examples
#' library(anomaly)
#' # generate some multivariate data
#' set.seed(0)
#' X<-simulate(n=1000,p=4,mu=10,locations=c(200,400,600),
#'             duration=100,proportions=c(0.25,0.5,0.75))
#' # compare the medians of each variate and transformed variate
#' head(apply(X,2,median))
#' head(apply(robustscale(X),2,median))
#' # compare the variances of each variate and transformed variate
#' head(apply(X,2,var))
#' head(apply(robustscale(X),2,var))
#'
#' @export
robustscale<-function(X)
{
    if(is.vector(X))
    {
        return((X-median(X))/mad(X))
    }
    else if(is.matrix(X))
    {
        return(Reduce(cbind,Map(function(i) array(robustscale(X[,i]),c(nrow(X),1)),1:ncol(X))))  
    }
    else
    {
        # incorrect type - throw an exception
    }
}
