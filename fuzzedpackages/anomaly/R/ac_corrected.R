#' Transforms the data X to account for autocorrelation.
#'
#' @name ac_corrected
#' 
#' @description Transforms the data X to account for autocorrelation by centring and scaling. It uses the transformation \eqn{ X_{i}^{'} = \frac{X_{i}-\mu_{i}}{k_{i}\sigma_{i}}}, were \eqn{\mu_{i}}
#' and \eqn{\sigma_{i}} are robust estimates for the mean and standard deviation of each variate (column), \eqn{X_{i}}, of X. The estimates are calculated using the median and median
#' absolute deviation. The scaling \eqn{k_{i} = \surd{\left( \frac{1+\phi_{i}}{1-\phi_{i}} \right)}}, with
#' \eqn{\phi_{i}} a robust estimate for the autocorrelation at lag 1, is used to account for AR(1) structure in the noise.
#' 
#' @param X A numeric matrix containing the potentially multivariate data to be transformed. Each column corresponds to a component and each row to an observation.
#' 
#' @return A numeric matrix of the same dimension as X containing the transformed data.
#'
#' @examples
#' library(anomaly)
#' # generate some multivariate data
#' set.seed(0)
#' X<-simulate(n=1000,p=4,mu=10,locations=c(200,400,600),
#'             duration=100,proportions=c(0.25,0.5,0.75))
#' # compare the medians of each variate and transformed variate
#' head(apply(X,2,median))
#' head(apply(ac_corrected(X),2,median))
#' # compare the variances of each variate and transformed variate
#' head(apply(X,2,var))
#' head(apply(ac_corrected(X),2,var))
#'
#' @export
ac_corrected<-function(X)
    {
        ac_corrections<-function(X)
        {
            X<-as.matrix(X)
            n<-dim(X)[1]
            m<-dim(X)[2]
            if(m == 1)
            {
	    	rcor<-covMcd(matrix(c(X[2:n],X[1:(n-1)]),ncol=2),cor=TRUE)
                psi<-rcor$cor[1,2]
                correction_factor<-sqrt((1-psi)/(1+psi))
                return(correction_factor)
            }
            else
            {
                return(unlist(Map(function(i) ac_corrections(as.matrix(X[,i])),1:m)))
            }
            
        }
        X<-robustscale(X)
        X<-t(ac_corrections(X)*t(X))
        return(X)
  }
  
