
.capa.mv.class<-setClass("capa.mv.class",contains="capa.class",representation())


capa.mv.class<-function(data,beta,beta_tilde,min_seg_len,max_seg_len,max_lag,type,
                     transform,anomaly_types,anomaly_positions,components,start_lags,end_lags,...)
{
.capa.mv.class(capa.class(data=data,beta=beta,beta_tilde=beta_tilde,min_seg_len=min_seg_len,max_seg_len=max_seg_len,max_lag=max_lag,type=type,
			transform=transform,anomaly_types=anomaly_types,anomaly_positions=anomaly_positions,components=components,start_lags=start_lags,end_lags=end_lags)
			,...)
}



#' @name point_anomalies
#'
#' @docType methods
#'
#' @rdname point_anomaly-methods
#'
#' @aliases point_anomalies,capa.mv.class-method
#'
#' @export
setMethod("point_anomalies",signature=list("capa.mv.class"),
          function(object)
          {
              return(callNextMethod(object))
          })


#' @name collective_anomalies
#'
#' @docType methods
#'
#' @rdname collective_anomalies-methods
#'
#' @aliases collective_anomalies,capa.mv.class-method
#'
#' @export
setMethod("collective_anomalies",signature=list("capa.mv.class"),
          function(object)
          {
              return(callNextMethod(object))
          })


#' @name summary
#'
#' @docType methods
#'
#' @rdname summary-methods
#'
#' @aliases summary,capa.mv.class-method
#' 
#' @export
setMethod("summary",signature=list("capa.mv.class"),function(object)
{
   summary(as(object,"capa.class"))
})


#' @name show
#'
#' @docType methods
#'
#' @rdname show-methods
#'
#' @aliases show,capa.mv.class-method
#' 
#' @export
setMethod("show",signature=list("capa.mv.class"),function(object)
{
  show(as(object,"capa.class"))
})


#' @name plot
#'
#' @docType methods
#'
#' @rdname plot-methods
#'
#' @aliases plot,capa.mv.class-method
#'
#' @export
setMethod("plot",signature=list("capa.mv.class"),function(x,subset,variate_names=FALSE,tile_plot)
{
    if(missing(subset))
    {
        subset<-1:ncol(x@data)
    }
    if(missing(tile_plot))
    {
        tile_plot<-NULL
    }
    return(plot(as(x,"capa.class"),subset=subset,variate_names=variate_names,tile_plot=tile_plot))
})



#'  Detection of multivariate anomalous segments and points using MVCAPA.
#'
#' @name capa.mv
#' 
#' @description This function implements MVCAPA (Multi-Variate Collective And Point Anomaly) from Fisch et al. (2019). 
#' It detects potentially lagged collective anomalies as well as point anomalies in multivariate time series data.  
#' The runtime of MVCAPA scales linearly (up to logarithmic factors) in \code{ncol(x)} and \code{maxlag}. If \code{max_seg_len} is not set, the runtime scales quadratically at worst and linearly 
#' at best in \code{nrow(x)}. If \code{max_seg_len} is set the runtime scales like \code{nrow(x)*max_seg_len}.
#' 
#' @param x A numeric matrix with n rows and p columns containing the data which is to be inspected.
#' @param beta A numeric vector of length p, giving the marginal penalties. If type ="meanvar" or if type = "mean"/"robustmean" and maxlag > 0 it defaults to the penalty regime 2' described in 
#' Fisch, Eckley, and Fearnhead (2019). If type = "mean"/"robustmean" and maxlag = 0 it defaults to the pointwise minimum of the penalty regimes 1, 2, and 3 in Fisch, Eckley, and Fearnhead (2019).
#' @param beta_tilde A numeric constant indicating the penalty for adding an additional point anomaly. It defaults to 3log(np), where n and p are the data dimensions.
#' @param type A string indicating which type of deviations from the baseline are considered. Can be "meanvar" for collective anomalies characterised by joint changes in mean and
#' variance (the default), "mean" for collective anomalies characterised by changes in mean only, or "robustmean" for collective anomalies characterised by changes in mean only which can be polluted by outliers.
#' @param min_seg_len An integer indicating the minimum length of epidemic changes. It must be at least 2 and defaults to 10.
#' @param max_seg_len An integer indicating the maximum length of epidemic changes. It must be at least the min_seg_len and defaults to Inf.
#' @param max_lag A non-negative integer indicating the maximum start or end lag. Default value is 0.
#' @param transform A function used to transform the data prior to analysis by \code{\link{capa.mv}}. This can, for example, be used to compensate for the effects of autocorrelation in the data. Importantly, the
#' untransformed data remains available for post processing results obtained using \code{\link{capa.mv}}. The package includes several methods that are commonly used for
#' the transform, (see \code{\link{robustscale}} and \code{\link{ac_corrected}}), but a user defined function can be specified. The default value is \code{transform=robust_scale}.
#' 
#' @return An instance of an S4 class of type capa.mv.class. 
#'
#' @references \insertRef{2019MVCAPA}{anomaly}
#'
#' @examples
#' library(anomaly)
#' 
#' ### generate some multivariate data
#' 
#' set.seed(0)
#' sim.data<-simulate(n=500,p=100,mu=2,locations=c(100,200,300),
#'                    duration=6,proportions=c(0.04,0.06,0.08))
#'                    
#' ### Apply MVCAPA
#' 
#' res<-capa.mv(sim.data,type="mean",min_seg_len=2)
#' plot(res)
#' 
#' ### generate some multivariate data
#' 
#' set.seed(2018)
#' x1 = rnorm(500)
#' x2 = rnorm(500)
#' x3 = rnorm(500)
#' x4 = rnorm(500)
#' 
#' ### Add two (lagged) collective anomalies
#' 
#' x1[151:200] = x1[151:200]+2
#' x2[171:200] = x2[171:200]+2
#' x3[161:190] = x3[161:190]-3
#' 
#' x1[351:390] = x1[371:390]+2
#' x3[351:400] = x3[351:400]-3
#' x4[371:400] = x4[371:400]+2
#' 
#' ### Add point anomalies
#'
#' x4[451] = x4[451]*max(1,abs(1/x4[451]))*5
#' x4[100] = x4[100]*max(1,abs(1/x4[100]))*5
#' x2[050] = x2[050]*max(1,abs(1/x2[050]))*5
#' 
#' my_x = cbind(x1,x2,x3,x4)
#' 
#' ### Now apply MVCAPA
#' 
#' res<-capa.mv(my_x,max_lag=20,type="mean")
#' 
#' plot(res)
#'
#' @export
capa.mv<-function(x,beta=NULL,beta_tilde=NULL,type="meanvar",min_seg_len=10,max_seg_len=Inf,max_lag=0,transform=robustscale)
{
    # data needs to be in the form of an array
    x<-to_array(x)
    if(dim(x)[2] < 2)
    {
      stop("data is not multivariate. Use capa or capa.uv for univariate analysis.")
    }
    res<-capa(x,beta,beta_tilde,type,min_seg_len,max_seg_len,max_lag,transform)
    return(
    capa.mv.class(data=res@data,
                 beta=res@beta,
	         beta_tilde=res@beta_tilde,
	         min_seg_len=res@min_seg_len,
	         max_seg_len=res@max_seg_len,
	         max_lag=as.integer(max_lag),
	         type=res@type,
                 transform=res@transform,
                 anomaly_types=res@anomaly_types,
	         anomaly_positions=res@anomaly_positions,
	         components=res@components,
	         start_lags=res@start_lags,
	         end_lags=res@end_lags)
           )
}


