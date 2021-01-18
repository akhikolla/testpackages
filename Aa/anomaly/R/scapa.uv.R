

.scapa.uv.class<-setClass("scapa.uv.class",contains="capa.class",representation())


scapa.uv.class<-function(data,beta,beta_tilde,min_seg_len,max_seg_len,max_lag,type,
                     transform,anomaly_types,anomaly_positions,components,start_lags,end_lags,...)
{
.scapa.uv.class(capa.class(data=data,beta=beta,beta_tilde=beta_tilde,min_seg_len=min_seg_len,max_seg_len=max_seg_len,max_lag=max_lag,type=type,
			transform=transform,anomaly_types=anomaly_types,anomaly_positions=anomaly_positions,components=components,start_lags=start_lags,end_lags=end_lags)
			,...)
}




#' @name point_anomalies
#'
#' @docType methods
#'
#' @rdname point_anomaly-methods
#'
#' @aliases point_anomalies,scapa.uv.class-method
#'
#' @export
setMethod("point_anomalies",signature=list("scapa.uv.class"),
          function(object,epoch=nrow(object@data))
          {
              return(callNextMethod(object,epoch=epoch)[,c(1,3)])
          })



#' @name collective_anomalies
#'
#' @docType methods
#'
#' @rdname collective_anomalies-methods
#'
#' @aliases collective_anomalies,scapa.uv.class-method
#'
#' @export
setMethod("collective_anomalies",signature=list("scapa.uv.class"),
          function(object,epoch=nrow(object@data))
          {
	      return(callNextMethod(object,epoch=epoch)[,c(1:2,6:7)])              
          })


#' @name plot
#'
#' @docType methods
#'
#' @param variate_name Logical value indicating if the variate name should be displayed. Default value is \code{variate_name=FALSE}.
#' 
#' @rdname plot-methods
#'
#' @aliases plot,scapa.uv.class-method
#'
#' @export
setMethod("plot",signature=list("scapa.uv.class"),function(x,epoch,variate_name=FALSE)
{
    if(missing(epoch))
    {
        epoch<-nrow(x@data)
    }
    return(plot(as(x,"capa.class"),epoch=epoch,variate_names=variate_name))
})





#' Detection of univariate anomalous segments using SCAPA.
#'
#' @name scapa.uv 
#'
#' @description An offline as-if-online implementation of SCAPA (Sequential Collective And Point Anomalies) by Bardwell et al. (2019) for online collective and point anomaly detection. This version of \code{capa.uv} has a default value
#' \code{transform=tierney} which uses sequential estimates for transforming the data prior to analysis. It also returns an S4 class which allows the results to be postprocessed
#' at different time points as if the data had been analysed in an online fashion up to that point.
#' 
#' @param x A numeric vector containing the data which is to be inspected.
#' @param beta A numeric vector of length 1 or \code{max_seg_len - min_seg_len + 1} indicating the penalty for adding additional collective anomalies of all possible
#' lengths. If an argument of length 1 is provided the same penalty is used for all collective anomalies irrespective of their length. The default value is 4log(n), where n denotes the number of observations.
#' @param beta_tilde A numeric constant indicating the penalty for adding an additional point anomaly. It defaults to 3log(n), where n is the number of observations.
#' @param type A string indicating which type of deviations from the baseline are considered. Can be "meanvar" for collective anomalies characterised by joint changes in mean and
#' variance (the default), "mean" for collective anomalies characterised by changes in mean only, or "robustmean" for collective anomalies characterised by changes in mean only which can be polluted by outliers.
#' @param min_seg_len An integer indicating the minimum length of epidemic changes. It must be at least 2 and defaults to 10.
#' @param max_seg_len An integer indicating the maximum length of epidemic changes. It must be at least the min_seg_len and defaults to Inf.
#' @param transform A function used to transform the data prior to analysis by \code{\link{scapa.uv}}. This can, for example, be used to compensate for the effects of autocorrelation in the data.
#' Importantly, the untransformed data remains available for post processing results obtained using \code{\link{scapa.uv}}. The package includes a method which can be used for
#' the transform, (see \code{\link{tierney}}, the default), but a user defined (ideally sequential) function can be specified.  
#'
#' @return An S4 class of type scapa.uv.class. 
#'
#' @references \insertRef{2018arXiv180601947F}{anomaly}
#' @references \insertRef{alex2020real}{anomaly}
#' 
#' @examples
#' library(anomaly)
#' 
#' # Simulated data example
#' # Generate data typically following a normal distribution with mean 0 and variance 1.
#' # Then introduce 3 anomaly windows and 4 point outliers.
#' 
#' set.seed(2018)
#' x  = rnorm(5000)
#' x[1601:1700] = rnorm(100,0,0.01)
#' x[3201:3300] = rnorm(100,0,10)
#' x[4501:4550] = rnorm(50,10,1)
#' x[c(1000,2000,3000,4000)] = rnorm(4,0,100)
#' # use magrittr to pipe the data to the transform
#' library(magrittr)
#' trans<-.%>%tierney(1000)
#' res<-scapa.uv(x,transform=trans)
#' 
#' # Plot results at two different times and note that anomalies are re-evaluated:
#' plot(res,epoch=3201)
#' plot(res,epoch=3205)
#' 
#' 
#' @export
scapa.uv<-function(x,beta=NULL,beta_tilde=NULL,type="meanvar",min_seg_len=10,max_seg_len=Inf,transform=tierney)
{
    # data needs to be in the form of an array
    x<-to_array(x)
    if(dim(x)[2] > 1)
    {
       stop("data for univariate analysis must have 1 variate. Use capa or capa.mv for multivariate data.")
    }
    res<-capa(x=x,beta=beta,beta_tilde=beta_tilde,type=type,min_seg_len=min_seg_len,max_seg_len=max_seg_len,transform=transform)
    return(
    scapa.uv.class(data=res@data,
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



