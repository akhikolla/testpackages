
.capa.uv.class<-setClass("capa.uv.class",contains="capa.class",representation())


capa.uv.class<-function(data,beta,beta_tilde,min_seg_len,max_seg_len,max_lag,type,
                     transform,anomaly_types,anomaly_positions,components,start_lags,end_lags,...)
{
.capa.uv.class(capa.class(data=data,beta=beta,beta_tilde=beta_tilde,min_seg_len=min_seg_len,max_seg_len=max_seg_len,max_lag=max_lag,type=type,
			transform=transform,anomaly_types=anomaly_types,anomaly_positions=anomaly_positions,components=components,start_lags=start_lags,end_lags=end_lags)
			,...)
}



#' @name point_anomalies
#'
#' @docType methods
#'
#' @rdname point_anomaly-methods
#'
#' @aliases point_anomalies,capa.uv.class-method
#'
#' @export
setMethod("point_anomalies",signature=list("capa.uv.class"),
          function(object)
          {
              return(callNextMethod(object)[,c(1,3)])
          })



#' @name collective_anomalies
#'
#' @docType methods
#'
#' @rdname collective_anomalies-methods
#'
#' @aliases collective_anomalies,capa.uv.class-method
#'
#' @export
setMethod("collective_anomalies",signature=list("capa.uv.class"),
          function(object)
          {
	      return(callNextMethod(object)[,c(1:2,6:7)])              
          })


#' @name show
#'
#' @docType methods
#'
#' @rdname show-methods
#'
#' @aliases show,capa.uv.class-method
#'
#' @export
setMethod("show",signature=list("capa.uv.class"),function(object)
{
  show(as(object,"capa.class"))
})

#' @name summary
#'
#' @docType methods
#'
#' @rdname summary-methods
#'
#' @aliases summary,capa.uv.class-method
#'
#' @export
setMethod("summary",signature=list("capa.uv.class"),function(object)
{
    summary(as(object,"capa.class"))
})


#' @name plot
#'
#' @docType methods
#'
#' @param variate_name Logical value indicating if the variate name should be displayed. Default value is \code{variate_name=FALSE}.
#' 
#' @rdname plot-methods
#'
#' @aliases plot,capa.uv.class-method
#'
#' @export
setMethod("plot",signature=list("capa.uv.class"),function(x,variate_name=FALSE)
{
    return(plot(as(x,"capa.class"),variate_names=variate_name))
})


# not exported - helper function for capa function
capa.uv_call<-function(x,beta=NULL,beta_tilde=NULL,type="meanvar",min_seg_len=10,max_seg_len=Inf)
{
    # configure defaults as required
    marshaller = marshall_MeanVarAnomaly
    if(type == "mean")
    {
        marshaller = marshall_MeanAnomaly
    }
    else if(type == "robustmean")
    {
      marshaller = marshall_RobustMeanAnomaly
    }
    if(is.null(beta))
    {
        if(type %in% c("mean","robustmean"))
        {
	    beta = 3*log(length(x))
        }
        else 
        {
	    beta = 4*log(length(x))
        }
    }
    if(length(beta) > 1 & length(beta) != (max_seg_len - min_seg_len + 1))
    {
        warning("beta has a number of entries larger than 1 but not equal to max_seg_len - min_seg_len + 1. Only the first one is kept.")
        beta = beta[1]
    }
    if(length(beta) == 1)
    {
        beta = rep(beta,max_seg_len - min_seg_len + 1)
    }
    if(is.null(beta_tilde))
    {
	beta_tilde = 3*log(length(x))
    }
    S<-marshaller(x,
                  as.integer(length(x)),
                  as.integer(min_seg_len),
                  as.integer(max_seg_len),
                  beta,
                  beta_tilde,
                  as.integer(1))
		     blob<-list(x,
               as.integer(length(x)),
               as.integer(min_seg_len),
               as.integer(max_seg_len),
               beta,
               beta_tilde,
               as.integer(1),
	       S)	       
    # construct the S4 capa class instance
    return(
	capa.class(array(x,c(length(x),1)), 
		     array(beta,c(length(beta),1)),
                     array(beta_tilde,c(1,1)),
                     as.integer(min_seg_len),
                     as.integer(max_seg_len),
                     integer(),
                     type,
                     function() return(),
                     S[seq(1,length(S),2)],
                     S[seq(2,length(S),2)],
                     array(1,c(length(x),1)), 
                     array(0,c(length(x),1)), 
                     array(0,c(length(x),1))) 
        )
}


#' Detection of univariate anomalous segments and points using CAPA.
#'
#' @name capa.uv
#' 
#' @description A technique for detecting anomalous segments and points in univariate time series data based on CAPA (Collective And Point Anomalies) by Fisch et al. (2018). CAPA assumes that the data has a certain mean and variance for most
#' time points and detects segments in which the mean and/or variance deviates from the typical mean and variance as collective anomalies. It also detects point
#' outliers and returns a measure of strength for the changes in mean and variance. If the number of anomalous windows scales linearly with the number of
#' data points, CAPA scales linearly with the number of data points. At
#' worst, if there are no anomalies at all and \code{max_seg_len} is unspecified, the computational cost of CAPA scales quadratically with the number of data points.
#'  
#' @param x A numeric vector containing the data which is to be inspected.
#' @param beta A numeric vector of length 1 or \code{max_seg_len - min_seg_len + 1} indicating the penalty for adding additional collective anomalies of all possible
#' lengths. If an argument of length 1 is provided the same penalty is used for all collective anomalies irrespective of their length. The default value is 4log(n), where n denotes the number of observations.
#' @param beta_tilde A numeric constant indicating the penalty for adding an additional point anomaly. It defaults to 3log(n), where n denotes the number of observations.
#' @param type A string indicating which type of deviations from the baseline are considered. Can be "meanvar" for collective anomalies characterised by joint changes in mean and
#' variance (the default), "mean" for collective anomalies characterised by changes in mean only, or "robustmean" for collective anomalies characterised by changes in mean only which can be polluted by outliers.
#' @param min_seg_len An integer indicating the minimum length of epidemic changes. It must be at least 2 and defaults to 10.
#' @param max_seg_len An integer indicating the maximum length of epidemic changes. It must be at least the min_seg_len and defaults to Inf.
#' @param transform A function used to transform the data prior to analysis by \code{\link{capa.uv}}. This can, for example, be used to compensate for the effects of autocorrelation
#' in the data. Importantly, the untransformed data remains available for post processing results obtained using \code{\link{capa.uv}}. The package includes several methods that are commonly used for
#' the transform, (see \code{\link{robustscale}} and \code{\link{ac_corrected}}), but a user defined function can be specified. The default values is \code{transform=robust_scale}. 
#'
#' @return An instance of an S4 class of type capa.uv.class. 
#'
#' @references \insertRef{2018arXiv180601947F}{anomaly}
#' 
#' @examples
#' library(anomaly)
#' data(machinetemp)
#' attach(machinetemp)
#' res<-capa.uv(temperature,type="mean")
#' canoms<-collective_anomalies(res)
#' dim(canoms)[1] # over fitted due to autocorrelation
#' psi<-0.98 # computed using covRob
#' inflated_penalty<-3*(1+psi)/(1-psi)*log(length(temperature))
#' res<-capa.uv(temperature,type="mean",beta=inflated_penalty,
#'              beta_tilde=inflated_penalty)
#' summary(res)
#' plot(res)
#'
#' library(anomaly)
#' data(Lightcurves)
#' ### Plot the data for Kepler 10965588: No transit apparent
#' plot(Lightcurves$Kepler10965588$Day,Lightcurves$Kepler10965588$Brightness,xlab = "Day",pch=".")
#' ### Examine a period of 62.9 days for Kepler 10965588
#' binned_data = period_average(Lightcurves$Kepler10965588,62.9)
#' inferred_anomalies = capa.uv(binned_data)
#' plot(inferred_anomalies)
#'
#' @export
capa.uv<-function(x,beta=NULL,beta_tilde=NULL,type="meanvar",min_seg_len=10,max_seg_len=Inf,transform=robustscale)
{
    # data needs to be in the form of an array
    x<-to_array(x)
    if(dim(x)[2] > 1)
    {
       stop("data for univariate analysis must have 1 variate. Use capa or capa.mv for multivariate data.")
    }
    res<-capa(x=x,beta=beta,beta_tilde=beta_tilde,type=type,min_seg_len=min_seg_len,max_seg_len=max_seg_len,transform=transform)
    return(
    capa.uv.class(data=res@data,
                 beta=res@beta,
	         beta_tilde=res@beta_tilde,
	         min_seg_len=res@min_seg_len,
	         max_seg_len=res@max_seg_len,
	         max_lag=res@max_lag,
	         type=res@type,
                 transform=res@transform,
                 anomaly_types=res@anomaly_types,
	         anomaly_positions=res@anomaly_positions,
	         components=res@components,
	         start_lags=res@start_lags,
	         end_lags=res@end_lags)
           )
}


