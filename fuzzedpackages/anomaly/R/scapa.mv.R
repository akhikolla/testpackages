

.scapa.mv.class<-setClass("scapa.mv.class",contains="capa.class",representation())


scapa.mv.class<-function(data,beta,beta_tilde,min_seg_len,max_seg_len,max_lag,type,
                     transform,anomaly_types,anomaly_positions,components,start_lags,end_lags,...)
{
.scapa.mv.class(capa.class(data=data,beta=beta,beta_tilde=beta_tilde,min_seg_len=min_seg_len,max_seg_len=max_seg_len,max_lag=max_lag,type=type,
			transform=transform,anomaly_types=anomaly_types,anomaly_positions=anomaly_positions,components=components,start_lags=start_lags,end_lags=end_lags)
			,...)
}



#' @name point_anomalies
#'
#' @docType methods
#'
#' @rdname point_anomaly-methods
#'
#' @aliases point_anomalies,scapa.mv.class-method
#'
#' @export
setMethod("point_anomalies",signature=list("scapa.mv.class"),
          function(object,epoch=nrow(object@data))
          {
              return(callNextMethod(object,epoch))
          })



#' @name collective_anomalies
#'
#' @docType methods
#'
#' @rdname collective_anomalies-methods
#'
#' @aliases collective_anomalies,scapa.mv.class-method
#'
#' @export
setMethod("collective_anomalies",signature=list("scapa.mv.class"),
          function(object,epoch=nrow(object@data))
          {
              return(callNextMethod(object,epoch=epoch))
          })


#' @name plot
#'
#' @docType methods
#'
#' @rdname plot-methods
#'
#' @aliases plot,scapa.mv.class-method
#'
#' @export
setMethod("plot",signature=list("scapa.mv.class"),function(x,subset,variate_names=FALSE,tile_plot,epoch)
{
    if(missing(epoch))
    {
        epoch<-nrow(x@data)
    }
    if(missing(subset))
    {
        subset<-1:ncol(x@data)
    }
    if(missing(tile_plot))
    {
        tile_plot<-NULL
    }
    return(plot(as(x,"capa.class"),subset=subset,variate_names=variate_names,tile_plot=tile_plot,epoch=epoch))
})



#' Online detection of multivariate anomalous segments and points using SMVCAPA.
#' 
#' @name scapa.mv 
#'
#' @description This function implements SMVCAPA from Fisch et al. (2019) in an as-if-online way. It detects potentially lagged collective anomalies as well as point anomalies in streaming data. 
#' The runtime scales linearly (up to logarithmic factors) in \code{ncol(x)}, \code{max_lag}, and \code{max_seg_len}. This version of \code{capa.uv} has a default value
#' \code{transform=tierney} which uses sequential estimates for transforming the data prior to analysis. It also returns an S4 class which allows the results to be postprocessed
#' as if the data had been analysed in an online fashion.
#' 
#' @param x A numeric matrix with n rows and p columns containing the data which is to be inspected.
#' @param beta A numeric vector of length p, giving the marginal penalties. If type ="meanvar" or if type = "mean" and maxlag > 0 it defaults to the penalty regime 2' described in 
#' Fisch, Eckley and Fearnhead (2019). If type = "mean" and maxlag = 0 it defaults to the pointwise minimum of the penalty regimes 1, 2, and 3 in Fisch, Eckley and Fearnhead (2019).
#' @param beta_tilde A numeric constant indicating the penalty for adding an additional point anomaly. It defaults to 3log(np), where n and p are the data dimensions. 
#' @param type A string indicating which type of deviations from the baseline are considered. Can be "meanvar" for collective anomalies characterised by joint changes in mean and
#' variance (the default), "mean" for collective anomalies characterised by changes in mean only, or "robustmean" for collective anomalies characterised by changes in mean only which can be polluted by outliers.
#' @param min_seg_len An integer indicating the minimum length of epidemic changes. It must be at least 2 and defaults to 10.
#' @param max_seg_len An integer indicating the maximum length of epidemic changes. It must be at least the min_seg_len and defaults to Inf.
#' @param max_lag A non-negative integer indicating the maximum start or end lag. Default value is 0.
#' @param transform A function used to transform the data prior to analysis by \code{\link{scapa.mv}}. This can, for example, be used to compensate for the effects of autocorrelation in the data. Importantly, the
#' untransformed data remains available for post processing results obtained using \code{\link{scapa.mv}}. The package includes a method which can be used for
#' the transform, (see \code{\link{tierney}}, the default), but a user defined (ideally sequential) function can be specified.  
#'
#' @return An S4 class of type scapa.mv.class. 
#' 
#' @references \insertRef{2019MVCAPA}{anomaly}
#' @references \insertRef{alex2020real}{anomaly}
#'
#' @examples
#' library(anomaly)
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
#' res<-scapa.mv(my_x,max_lag=20,type="mean")
#' 
#' ### Examine the output at different times and see how the results are updated:
#' 
#' plot(res,epoch=155)
#' plot(res,epoch=170)
#' plot(res,epoch=210)
#'
#' @export
scapa.mv<-function(x,beta=NULL,beta_tilde=NULL,type="meanvar",min_seg_len=10,max_seg_len=Inf,max_lag=0,transform=tierney)
{
    # data needs to be in the form of an array
    x<-to_array(x)
    if(dim(x)[2] < 2)
    {
      stop("data is not multivariate. Use capa or capa.uv for univariate analysis.")
    }
    res<-capa(x,beta,beta_tilde,type,min_seg_len,max_seg_len,max_lag,transform)
    return(
    scapa.mv.class(data=res@data,
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