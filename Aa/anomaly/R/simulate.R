#' A function for generating simulated multivariate data
#'
#' @name simulate
#'
#' @description Generates multivariate simulated data having n observations and p variates. The data have a standard Gaussian distribution except at
#' a specified number of locations where there is a change in mean in a proportion of the variates. The function is useful for generating data to demonstrate and assess
#' multivariate anomaly detection methods such as \code{capa.mv} and \code{pass}.
#' 
#' @param n The number of observations. The default is \code{n=100}.
#' @param p The number of variates. The default is \code{p=10}.
#' @param mu The change in mean. Default is \code{mu=1}.
#' @param locations A vector of locations (or scalar for a single location) where the change in mean occurs. The default is \code{locations=20}.
#' @param durations A scalar or vector (the same length as \code{locations}) of values indicating the duration for the change in mean. If the durations are all
#' of the same length then a scalar value can be used. The default is \code{durations=20}.
#' @param proportions A scalar or vector (the same length as \code{locations}) of values in the range (0,1] indicating the proportion of variates at each location that are affected by
#' the change in mean. If the proportions are all same than a scalar value can be used. The default is \code{proportions=0.1}.
#' 
#' @return A matrix with n rows and p columns
#'
#'
#' @examples
#' library(anomaly)
#' sim.data<-simulate(500,200,2,c(100,200,300),6,c(0.04,0.06,0.08))
#' 
#' @export
simulate<-function(n=100,p=10,mu=1,locations=40,durations=20,proportions=0.1)
{
    if(length(n) > 1)
    {
        stop("n must be a scalar")
    }
    if(length(p) > 1)
    {
        stop("p must be a scalar")
    }
    if(length(mu) > 1)
    {
        stop("mu must be a scalar")
    }
    if(n < 1)
    {
        stop("n must be > 0")
    }
    if(p < 1)
    {
        stop("p must be > 0")
    }
    if(!Reduce("&&",locations > 0))
    {
        stop("location values must all be > 0")   
    }
    if(!Reduce("&&",durations > 0))
    {
        stop("durations must all be > 0")   
    }
    if(!Reduce("&&",proportions > 0))
    {
        stop("proportion values  must all be > 0")   
    }
    if(!Reduce("&&",proportions <= 1))
    {
        stop("proportion values  must all be > 0")
    }
    if(!Reduce("&&",proportions <= 1))
    {
        stop("proportion values  must all be > 0")
    }
    if(length(proportions) != 1 && length(proportions) != length(locations))
    {
        stop("proportions must be a scalar or a vector the same size as locations")
    }
    if(length(durations) != 1 && length(durations) != length(locations))
    {
        stop("durations must be a scalar or a vector the same size as locations")
    }
    if(length(durations) == 1)
    {
        durations<-rep(durations,length(locations))
    }
    if(length(proportions) == 1)
    {
        proportions<-rep(proportions,length(locations))
    }
    if(!Reduce("&",locations+durations < n))
    {
        stop("locations+durations must be < n")
    }
    
    q = length(proportions);
    if(length(durations) == 1)
        {
            s = rep(durations, q);
        }
    else
        {
            s=durations
        }
    X = matrix(0, p, n);
    for (i in 1 : p)
    {
        X[i,] = rnorm(n);
    }
    for (j in 1:q)
    {
      # check if proportion of series affected is non zero
      paffected = round(proportions[j]*p)
      if (paffected > 0){
        for (i in 1:paffected)
        {
          X[i,locations[j]:(locations[j]+s[j]-1)] = X[i,locations[j]:(locations[j]+s[j]-1)]+mu;
        }
      }
    }
    return(t(X))
}
