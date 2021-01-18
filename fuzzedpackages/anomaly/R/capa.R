.capa.class<-setClass("capa.class",representation(data="array",beta="array",beta_tilde="vector",min_seg_len="integer",max_seg_len="integer",max_lag="integer",type="character",
                                                  transform="function",anomaly_types="vector",anomaly_positions="vector",components="array",start_lags="array",end_lags="array"))

capa.class<-function(data,beta,beta_tilde,min_seg_len,max_seg_len,max_lag,type,
                     transform,anomaly_types,anomaly_positions,components,start_lags,end_lags,...)
{
    .capa.class(data=data,beta=beta,beta_tilde=beta_tilde,min_seg_len=min_seg_len,max_seg_len=max_seg_len,max_lag=max_lag,type=type,
                transform=transform,anomaly_types=anomaly_types,anomaly_positions=anomaly_positions,components=components,start_lags=start_lags,end_lags=end_lags,...)
}


# utility function to coerce data to an array structure
to_array<-function(X)
{
  if(is.data.frame(X))
  {
     X<-data.matrix(X, rownames.force = NA)
  }
  X<-as.array(X)
  dims<-dim(X)
  if(length(dims) == 1)
  {
    X<-array(X,c(dims,1))
  }
  if(length(dims) > 2)
  {
    stop("data in array structures with dimension > 2 not supported") 
  }
  return(X)
}

# utility function to process commen features of summary and show
summary_show_common<-function(object,epoch=nrow(object@data))
{
    if(epoch < 0)
    {
        stop("epoch should be a positive integer")
    }
    if(epoch > nrow(object@data))
    {
        stop("epoch cannot be greater than the number of observations in the data")
    }
    if(dim(object@data)[2] == 1)
    {
       cat("Univariate ",sep="")
    }
    else
    {
       cat("Multivariate ",sep="")	
    }
    cat("CAPA detecting changes in ",sep="")
    if(object@type == "meanvar")
    {
        cat("mean and variance.","\n",sep="") 
    }
    if(object@type %in% c("mean","robustmean"))
    {
        cat("mean.","\n",sep="") 
    }
    cat("observations = ",dim(object@data)[1],sep="")
    cat("\n",sep="")
    if(dim(object@data)[2] != 1)
    {
       cat("variates = ",dim(object@data)[2],"\n",sep="")
    }
    cat("minimum segment length = ",object@min_seg_len,'\n',sep="")
    cat("maximum segment length = ",object@max_seg_len,'\n',sep="")
    if(dim(object@data)[2] != 1)
    {
       cat("maximum lag = ",object@max_lag[1],'\n',sep="")
    }
    if(epoch != nrow(object@data))
    {
       cat("epoch = ",epoch,"\n",sep="")
    }
}


#' Point anomaly location and strength.
#'
#' @name point_anomalies
#'
#' @description Creates a data frame containing point anomaly locations and strengths as detected by \code{\link{capa}}, \code{\link{capa.uv}}, \code{\link{capa.mv}}, \code{\link{scapa.uv}}, and \code{\link{scapa.mv}}.
#' 
#'
#' For an object produced by \code{\link{capa.uv}} or \code{\link{scapa.uv}}, the output is a data frame  with columns containing the position and
#' strength of the anomaly. 
#' 
#' For an object produced by \code{\link{capa.mv}} or \code{\link{capa.mv}}, \code{point_anomalies} returns a data frame with columns containing the position, variate, and
#' strength of the anomaly. 
#'
#' 
#' For an object produced by \code{\link{capa}}, \code{point_anomalies} returns the same results as \code{\link{scapa.uv}} when the data is univariate, and the same results as
#' \code{\link{scapa.mv}} when the data is multivariate.
#'
#' 
#' @docType methods
#'
#' @rdname point_anomaly-methods
#'
#' @seealso \code{\link{capa}}, \code{\link{capa.uv}}, \code{\link{capa.mv}}.
#'
if(!isGeneric("point_anomalies")) {setGeneric("point_anomalies",function(object,...) {standardGeneric("point_anomalies")})}

#' @name point_anomalies
#' @param object An instance of an S4 class produced by \code{\link{capa}}, \code{\link{capa.uv}}, and \code{\link{capa.mv}}.
#' @param epoch Positive integer. CAPA methods are sequential and as such, can generate results up to, and including, any epoch within the data series. This can be controlled by the value
#' of \code{epoch} and is useful for examining how the inferred anomalies are modified as the data series grows. The default value for \code{epoch} is the length of the data series.
#' 
#' @return A data frame. 
#'
#' @rdname point_anomaly-methods
#'
#' @aliases point_anomalies,capa.class-method
#'
#' @export
setMethod("point_anomalies",signature=list("capa.class"),
          function(object,epoch=nrow(object@data))
          {
              if(epoch < 0)
              {
                  stop("epoch should be a positive integer")
              }
              if(epoch > nrow(object@data))
              {
                  stop("epoch cannot be greater than the number of observations in the data")
              }
              # get the anomalies
              anoms<-anomalies(object,epoch)
              # transform data
              data_dash<-object@transform(object@data)
              p_anoms<-Map(function(x) x[1],Filter(function(x) x[1] == x[2],anoms))
              if(length(p_anoms) > 0)
              {
                  p_anom_daf = Reduce(rbind,
                       Map(
                         function(p_anom)
                         {
                           variates<-seq(1,ncol(data_dash))[as.logical(object@components[p_anom[1],])]
                           location<-rep(p_anom[1],length(variates))
                           strength<-abs(data_dash[p_anom[1],variates])
                           return(
                             data.frame("location"=location,
                                        "variate"=variates,
                                        "strength"=strength)
                           )
                         },
                         p_anoms)
                  )
              }
              
              extra_anoms = data.frame("location"=integer(0),"variate"=integer(0),"strength"=integer(0))
              
              if (object@type == "robustmean"){
                
                tmp = collective_anomalies(as(object,"capa.class"))
                
                if (nrow(tmp)>0)
                {
                  
                  extra_anoms = Reduce(rbind,
                                      Map(
                                        function(ii)
                                        {
                                          relevant_row = tmp[ii,]
                                          if (is.null(relevant_row$start.lag)){
                                            effective_start = relevant_row$start
                                            effective_end   = relevant_row$end
                                          } else{
                                            effective_start = relevant_row$start+relevant_row$start.lag
                                            effective_end   = relevant_row$end-relevant_row$start.lag                                      
                                          }
                                          x_data = data_dash[effective_start:effective_end,relevant_row$variate]
                                          standardised_x_data = x_data - tukey_mean(x_data,sqrt(object@beta_tilde))
                                          
                                          location<-which(abs(standardised_x_data)>sqrt(object@beta_tilde))
                                          strength<-abs(standardised_x_data[location])
                                          variates<-rep(relevant_row$variate,length(location))
                                          if (length(location > 0)){
                                            location = location - 1 + effective_start
                                          }
                                          return(
                                            data.frame("location"=location,
                                                       "variate"=variates,
                                                       "strength"=strength)
                                          )
                                        },
                                        1:nrow(tmp))
                  )
                  
                }
              
              }
              
              if(length(p_anoms) + nrow(extra_anoms) == 0)
              {
                  return(data.frame("location"=integer(0),"variate"=integer(0),"strength"=integer(0)))
              }
              else
              {
                  out = rbind(p_anom_daf,extra_anoms)
                  return(out[order(out$location,out$variate),])
              }
          }
          )




# helper for collective_anomaies with object of type capa.class
merge_collective_anomalies<-function(object,epoch)
{
   if(object@type != "mean")
    {
        stop("merge=TRUE option only available when type=mean")
    }
   canoms<-collective_anomalies(object,merged=FALSE)
   # get the location (start,end) for each collective anomaly
   locations<-unique(canoms[,1:2])
   # get the sum of the test statistic over the variates affected by each collective anomaly
   sums<-unlist(Map(function(start) sum(canoms[canoms[,1] == start,]$mean.change),locations[,1]))
   # order the locations
   locations<-locations[order(sums,decreasing=TRUE),]
   # sort the sums
   sums<-sort(sums,decreasing=TRUE)
   # divide sums by the sum of the betas
   sums<-sums/sum(object@beta)
   return(data.frame("start"=locations[,1],"end"=locations[,2],"mean.change"=sums))
}



#' Collective anomaly location, lags, and mean/variance changes.
#'
#' @name collective_anomalies
#'
#' @description Creates a data frame containing collective anomaly locations, lags and changes in mean and variance as detected by
#' \code{\link{capa.uv}}, \code{\link{capa.mv}}, \code{\link{scapa.uv}}, \code{\link{scapa.mv}}, \code{\link{capa}}, \code{\link{pass}}, and \code{\link{sampler}}. 
#'
#' For an object produced by \code{\link{capa.uv}} or \code{\link{scapa.uv}}, \code{collective_anomalies} returns a data frame with columns containing the start and end position of the anomaly, the change in mean
#' due to the anomaly. When \code{type="meanvar"}, the change in variance due to the anomaly is also returned in an additional column.
#' 
#'
#' For an object produced by \code{\link{capa.mv}} or \code{\link{scapa.mv}}, \code{collective_anomalies} returns a data frame with columns containing the start and end position of the anomaly, the variates 
#' affected by the anomaly, as well as their the start and end lags. When \code{type="mean"/"robustmean"} only the change in mean is reported. When \code{type="meanvar"} both the change in mean and
#' change in variance are included. If \code{merged=FALSE} (the default), then all the collective anomalies are processed individually even if they are common across multiple variates.
#' If \code{merged=TRUE}, then the collective anomalies are grouped together across all variates that they appear in.
#'
#'
#' For an object produced by \code{\link{capa}}, \code{collective_anomalies} returns the same results as \code{\link{scapa.uv}} when the data is univariate, or the same results as
#' \code{\link{scapa.mv}} when the data is multivariate.
#'
#' For an object produced by \code{\link{pass}} or \code{sampler} returns a data frame containing the start, end and strength of the collective anomalies.
#'
#' 
#' @docType methods
#'
#' @rdname collective_anomalies-methods
#'
if(!isGeneric("collective_anomalies")) {setGeneric("collective_anomalies",function(object,...) {standardGeneric("collective_anomalies")})}

#' @name collective_anomalies
#' @param object An instance of an S4 class produced by \code{\link{capa}}, \code{\link{capa.uv}} and \code{\link{capa.mv}}.
#' @param epoch Positive integer. CAPA methods are sequential and as such, can generate results up to, and including, any epoch within the data series. This can be controlled by the value
#' of \code{epoch} and is useful for examining how the inferred anomalies are modified as the data series grows. The default value for \code{epoch} is the length of the data series.
#' @param merged Boolean value. If \code{merged=TRUE} then collective anomalies that are common across multiple variates are merged together. This is useful when comparing the relative strength
#' of multivariate collective anomalies. Default value is \code{merged=FALSE}. Note - \code{merged=TRUE} is currently only available when \code{type="mean"}.  
#' 
#' @return A data frame.
#' 
#' @rdname collective_anomalies-methods
#'
#' @aliases collective_anomalies,capa.class-method
#' 
#' @seealso \code{\link{capa}},\code{\link{capa.uv}},\code{\link{capa.mv}}. 
#'
#' @export
setMethod("collective_anomalies",signature=list("capa.class"),
          function(object,epoch=nrow(object@data),merged=FALSE)
          {
              if(epoch < 0)
              {
                  stop("epoch should be a positive integer")
              }
              if(epoch > nrow(object@data))
              {
                  stop("epoch cannot be greater than the number of observations in the data")
              }
              if(merged)
              {
                  return(merge_collective_anomalies(object,epoch))
              }
              # get the anomalies
              anoms<-anomalies(object,epoch)
              # transform data
              data_dash<-object@transform(object@data)
              c_anoms<-Filter(function(x) x[1] != x[2],anoms)
              if(length(c_anoms) == 0)
              {                      
                  return(data.frame("start"=integer(0),"end"=integer(0),"variate"=integer(0),"start.lag"=integer(0),"end.lag"=integer(0),"mean.change"=integer(0),"variance.change" = integer(0)))
              }
              res<-Reduce(rbind,
                          Map(
                              function(c_anom)
                              {
                                  variates<-seq(1,ncol(data_dash))[as.logical(object@components[c_anom[2],])]
                                  start<-rep(c_anom[1],length(variates))
                                  end<-rep(c_anom[2],length(variates))
                                  start_lags<-object@start_lags[c_anom[2],variates]
                                  end_lags<-object@end_lags[c_anom[2],variates]
                                  return(
                                      data.frame("start"=start,
                                                 "end"=end,
                                                 "variate"=variates,
                                                 "start.lag"=start_lags,
                                                 "end.lag"=end_lags)
                                  )
                              },
                              c_anoms)
                          )
              if(object@type == "meanvar")
              {
                  changes<-data.frame(
                      Reduce(rbind,
                             Map(
                                 function(variate,start,end,start_lag,end_lag)
                                 {
                                     variance<-var(data_dash[(start+start_lag):(end-end_lag),variate])
                                     variance_change<-sqrt(variance)+1/sqrt(variance)-2
                                     mean_change<-mean(data_dash[(start+start_lag):(end-end_lag),variate])^2/sqrt(variance)
                                     return(array(c(mean_change,variance_change),c(1,2)))    
                                 },
                                 res$variate,
                                 res$start,
                                 res$end,
                                 res$start.lag,
                                 res$end.lag)
                             ),row.names=NULL
                  )
                  names(changes)<-c("mean.change","variance.change")
              }
              else if(object@type %in% c("mean","robustmean"))
              {
                  changes<-data.frame(
                      Reduce(rbind,
                             Map(
                                 function(variate,start,end,start_lag,end_lag)
                                 {
                                     if (object@type == "mean"){
                                          mean_change<-mean(data_dash[(start+start_lag):(end-end_lag),variate])^2
                                     }
                                     if (object@type == "robustmean"){
                                          mean_change<-tukey_mean(data_dash[(start+start_lag):(end-end_lag),variate],sqrt(object@beta_tilde))^2
                                     }
                                     variance_change<-mean_change*((end-end_lag)-(start+start_lag)+1) 
                                     return(array(c(mean_change,variance_change),c(1,2)))
                                     # return(array(c(mean_change),c(1,1)))    
                                 },
                                 res$variate,
                                 res$start,
                                 res$end,
                                 res$start.lag,
                                 res$end.lag)
                             ),row.names=NULL
                  )
                  names(changes)<-c("mean.change","test.statistic")
                  # names(changes)<-c("mean.change")
              }
              else
              {
                  # default - no summary information produced
                  return(res)
              }
              res<-cbind(res,changes);
              # tighten
	      start.lag<-end.lag<-"."<-"%<>%"<-NULL # circumvent CRAN check issues
	      res %<>% group_by(.,start) %>% mutate(.,end.lag=end.lag-min(end.lag),start.lag=start.lag-min(start.lag)) %>% as.data.frame(.)
              return(res)
          }
          )







#' Summary of collective and point anomalies.
#'
#' @name summary
#'
#' @description Summary methods for S4 objects returned by \code{\link{capa}}, \code{\link{capa.uv}}, \code{\link{capa.mv}},
#' \code{\link{pass}}, and \code{\link{sampler}}.  The output displayed
#' depends on the type of object passed to summary. For all types, the output indicates whether the data is univariate or
#' multivariate, the number of observations in the data, and the type of change being detected.
#'
#'
#' For an object produced by \code{\link{capa.uv}} or \code{\link{capa.mv}}, \code{\link{pass}}, or \code{\link{sampler}},  \code{summary}
#' displays a summary of the analysis.
#' 
#' For an object produced by \code{\link{capa}} is the same as for an object produced by \code{\link{capa.uv}}
#' or \code{\link{capa.mv}}. 
#'
#' @docType methods
#'
#' @param object An instance of an S4 class produced by \code{\link{capa}}, \code{\link{capa.uv}}, \code{\link{capa.mv}}, or \code{\link{pass}}.
#' @param epoch Positive integer. CAPA methods are sequential and as such, can generate results up to, and including, any epoch within the data series. This can be controlled by the value
#' of \code{epoch} and is useful for examining how the inferred anomalies are modified as the data series grows. The default value for \code{epoch} is the length of the data series.
#' @param ... Ignored.
#'
#' @rdname summary-methods
#'
#' @aliases summary,capa.class-method
#' 
#' @seealso \code{\link{capa}},\code{\link{capa.uv}},\code{\link{capa.mv}},\code{\link{pass}},\code{\link{sampler}}. 
#'
#' @export
setMethod("summary",signature=list("capa.class"),function(object,epoch=nrow(object@data))
{
  summary_show_common(object,epoch)
  p_anoms<-point_anomalies(object,epoch)
  c_anoms<-collective_anomalies(object,epoch)
  cat("\n",sep="")
  cat("Point anomalies detected : ",nrow(p_anoms),"\n",sep="")
  if (nrow(p_anoms)>0)
  {
     print(p_anoms)
  }
  cat("\n",sep="")
  cat("Collective anomalies detected : ",length(unique(c_anoms$start)),"\n",sep="")
  if (nrow(c_anoms)>0)
  {
    print(c_anoms)
  }
  invisible()
})



#' Displays S4 objects produced by capa methods.
#'
#' @name show
#'
#' @description Displays S4 object produced by \code{\link{capa}}, \code{\link{capa.uv}}, \code{\link{capa.mv}}, \code{\link{pass}}, \code{\link{bard}}, and \code{\link{sampler}}.
#' The output displayed depends on the type of S4 object passed to the method. For all types, the output indicates whether the data is univariate or
#' multivariate, the number of observations in the data, and the type of change being detected.
#'
#' @docType methods
#'
#' @param object An instance of an S4 class produced by \code{\link{capa}}, \code{\link{capa.uv}}, \code{\link{capa.mv}}, \code{\link{pass}}, \code{\link{bard}}, or \code{\link{sampler}}.
#'
#' @rdname show-methods
#'
#' @aliases show,capa.class-method
#' 
#' @seealso \code{\link{capa}},\code{\link{capa.uv}},\code{\link{capa.mv}},\code{\link{pass}},\code{\link{bard}},\code{\link{sampler}}. 
#'
#' @export
setMethod("show",signature=list("capa.class"),function(object)
{
    summary_show_common(object)
    p_anoms<-point_anomalies(object)
    c_anoms<-collective_anomalies(object)
    cat("Point anomalies detected : ",nrow(p_anoms),"\n",sep="")
    cat("Collective anomalies detected : ",length(unique(c_anoms$start)),"\n",sep="")
    invisible()
})


anomalies<-function(x,epoch=NULL)
{
    anomaly_types<-x@anomaly_types
    anomaly_positions<-x@anomaly_positions
    if(is.null(epoch))
    {
        epoch<-length(anomaly_types) 
    }
    if(epoch < 0 || epoch > length(anomaly_types))
    {
	stop("epoch must be a positive integer <= number of observations in data")
    }
    anomaly_types<-anomaly_types[1:epoch]
    anomaly_positions<-anomaly_positions[1:epoch]
    anom<-list()
    cursor<-length(anomaly_positions)
    count<-0
    while(cursor > 1)
    {
        if(anomaly_types[cursor] == 1)
        {
            count<-count+1
            anom[[count]]<-c(cursor,cursor)
        }
        if(anomaly_types[cursor] == 2)
        {
            count<-count+1
            anom[[count]]<-c(anomaly_positions[cursor]+1,cursor)
        }
        cursor<-anomaly_positions[cursor]
    }
    return(rev(anom))
}


 

# not exported - helper function used by capa function
capa.mv_call<-function(x,beta=NULL,beta_tilde=NULL,type="meanvar",min_seg_len=10,max_seg_len=Inf,max_lag=0)
{
    # configure defaults as required
    marshaller = marshall_MeanVarAnomalyMV
    if(type == "mean")
    {
        marshaller = marshall_MeanAnomalyMV
    }
    if(type == "robustmean")
    {
        marshaller = marshall_RobustMeanAnomalyMV
    }
    # process beta
    n = nrow(x)
    p = ncol(x)
    s = 1.5*log(n)
    if (is.null(beta))
    {            
        if(type %in% c("mean","robustmean"))
        {
            if (max_lag == 0)
            {
                to_be_tried = 1:p
                a_vector = qchisq(seq(p-1,0)/p, 1)
        
                penalty_1 = 2*s + 2*to_be_tried*log(p)
                penalty_2 = rep(p + 2*s + 2*sqrt(p*s),p)
        
                penalty_3    = 2*(s+log(p)) + 1:p + 2*p*a_vector*dchisq(a_vector, 1) + 2*sqrt((1:p+2*p*a_vector*dchisq(a_vector, 1))*(s+log(p)))
                penalty_3[p] = 2*(s+log(p)) + p + 2*sqrt(p*(s+log(p)))
        
                beta = diff( c(0,pmin(penalty_1,penalty_2,penalty_3)))   
            }
      
            if (max_lag > 0)
            {
                beta = rep(2*log(p*(max_lag+1)),p)
                beta[1] = beta[1] + 2*s
            }
        }
        if(type == "meanvar")
        {
            beta = rep(4*log(p*(max_lag+1)),p)
            beta[1] = beta[1] + 4*s         
        }
    }
    else if(length(beta) == 1)
    {
        beta = rep(beta,p)
    }
    if(is.null(beta_tilde))
    {
        beta_tilde = 3*log(length(x))
    }
    
    # call the underlying method
    S<-marshaller(as.vector(x),
                  as.integer(nrow(x)),
                  as.integer(ncol(x)),
                  as.integer(max_lag),
                  as.integer(min_seg_len),
                  beta,
                  beta_tilde,
                  as.integer(max_seg_len),
                  as.integer(1))
    # construct the S4 capa class instance
    S = matrix(S,nrow(x),byrow=T)
    p = ncol(x)
    return(
        capa.class(x,
                   array(beta,c(length(beta),1)),
                   array(beta_tilde,c(1,1)),
                   as.integer(min_seg_len),
                   as.integer(max_seg_len),
                   integer(max_lag),
                   type,
                   function() return(),
                   S[,1], 
                   S[,2],
                   matrix(S[,2 + 0*p + 1:p],ncol=p),
                   matrix(S[,2 + 1*p + 1:p],ncol=p),  
                   matrix(S[,2 + 2*p + 1:p],ncol=p))
        )
}

#' A technique for detecting anomalous segments and points based on CAPA.
#'
#' @name capa 
#'
#' @description A technique for detecting anomalous segments and points based on CAPA (Collective And Point Anomalies) by Fisch et al. (2018). This is a generic method that can be used for both univariate
#' and multivariate data. The specific method that is used for the analysis is deduced by \code{capa} from the dimensions of the data.
#' 
#' @param x A numeric matrix with n rows and p columns containing the data which is to be inspected.
#' @param beta A numeric vector of length p, giving the marginal penalties. If p > 1, type ="meanvar" or type = "mean" and max_lag > 0 it defaults to the penalty regime 2' described in 
#' Fisch, Eckley and Fearnhead (2019). If p > 1, type = "mean"/"meanvar" and max_lag = 0 it defaults to the pointwise minimum of the penalty regimes 1, 2, and 3 in Fisch, Eckley and Fearnhead (2019).
#' @param beta_tilde A numeric constant indicating the penalty for adding an additional point anomaly. It defaults to 3log(np), where n and p are the data dimensions.
#' @param type A string indicating which type of deviations from the baseline are considered. Can be "meanvar" for collective anomalies characterised by joint changes in mean and
#' variance (the default), "mean" for collective anomalies characterised by changes in mean only, or "robustmean" for collective anomalies characterised by changes in mean only which can be polluted by outliers.
#' @param min_seg_len An integer indicating the minimum length of epidemic changes. It must be at least 2 and defaults to 10.
#' @param max_seg_len An integer indicating the maximum length of epidemic changes. It must be at least min_seg_len and defaults to Inf.
#' @param max_lag A non-negative integer indicating the maximum start or end lag. Only useful for multivariate data. Default value is 0.
#' @param transform A function used to centre the data prior to analysis by \code{\link{capa}}. This can, for example, be used to compensate for the effects of autocorrelation in the data.
#' Importantly, the untransformed data remains available for post processing results obtained using \code{\link{capa}}. The package includes several methods that are commonly used for
#' the transform, (see \code{\link{robustscale}} and \code{\link{ac_corrected}}), but a user defined function can be specified. The default values is \code{transform=robust_scale}. 
#' 
#' @return An instance of an S4 class of type capa.class. 
#'
#' @references \insertRef{2018arXiv180601947F}{anomaly}
#'
#'
#' @examples
#' library(anomaly)
#' # generate some multivariate data
#' set.seed(0)
#' sim.data<-simulate(n=500,p=100,mu=2,locations=c(100,200,300),
#'                    duration=6,proportions=c(0.04,0.06,0.08))
#' res<-capa(sim.data,type="mean",min_seg_len=2,max_lag=5)
#' collective_anomalies(res)
#' plot(res)
#' 
#' @export
#'
capa<-function(x,beta=NULL,beta_tilde=NULL,type="meanvar",min_seg_len=10,max_seg_len=Inf,max_lag=0,transform=robustscale)
{
    # make sure the callable transform object is of type function
    # data needs to be in the form of an array
    x<-to_array(x)
    if(!is_array(x))
    {
        stop("cannot convert x to an array")
    }
    if(!all(is_not_na(x)))
    {
        stop("x contains NA values")
    }
    if(!all(is_not_null(x)))
    {
        stop("x contains NULL values")
    }
    if(!is_numeric(x))
    {
        stop("x must be of type numeric")
    }
    if(Reduce("|",x == Inf))
    {
        stop("x contains Inf values")
    }
    # check dimensions
    if(length(dim(x)) != 2)
    {
        stop("cannot convert x to a two dimensional array")
    }
    if(dim(x)[1] == 1)
    {
        x<-t(x)
    }
    # try transforming the data
    if(!is_function(transform))
    {
        stop("transform must be a function")
    }
    x_untransformed<-x
    x<-transform(x)
    # and check the transformed data
    if(!is_array(x))
    {
        stop("cannot convert the transformed x to an array - check the transform function")
    }
    if(!all(is_not_na(x)))
    {
        stop("the transformed x contains NA values - check the transform function")
    }
    if(!all(is_not_null(x)))
    {
        stop("the transformed x contains NULL values - check the transform function")
    }
    if(!is_numeric(x))
    {
        stop("the transformed x must be of type numeric - check the transform function")
    }
    if(length(dim(x)) != 2)
    {
        stop("cannot convert transformed x to a two dimensional array - check the transform function")
    }

    # check the type 
    types=list("mean","robustmean","meanvar")
    if(!(type %in% types))
    {
        stop("type has to be one of mean, robustmean, or meanvar")
    }

    # set analysis type
    univariate<-FALSE
    if(dim(x)[2] == 1)
    {
        univariate<-TRUE
    }

    # check max_lag
    if(!is_numeric(max_lag))
    {
        stop("max_lag must be numeric")
    }
    if(max_lag < 0)
    {
        stop("max_lag must be positive or zero")
    }
    if(univariate && max_lag != 0)
    {
        warning("max_lag != 0 redundant for a univariate analysis")
    }

    # check min_seg_len
    if(!is_numeric(min_seg_len))
    {
        stop("min_seg_len must be numeric")
    }
    if(type %in% c("mean","robustmean") && min_seg_len < 1)
    {
        stop("when type=mean or type=robustmean, min_seg_len must be greater than zero")
    }
    if(type=="meanvar" && min_seg_len < 2)
    {
        stop("when type=meanvar, min_seg_len must be greater than 1")
    }
    if(min_seg_len < 1)
    {
        stop("min_seg_len must be greater than zero")
    }
    if(min_seg_len > dim(x)[1])
    {
        stop("min_seg_len must be less tha the number of observations in x")
    }

    # check max_seg_len
    if(is.null(max_seg_len))
    {
        max_seg_len=dim(x)[1]
    }
    if(max_seg_len == Inf)
    {
        max_seg_len = length(x)
    }
    if(!is_numeric(max_seg_len))
    {
        stop("max_seg_len must be numeric")
    }
    if(max_seg_len < 1)
    {
        stop("max_seg_len must be greater than zero")
    }

    # check relative values of min and max segment length
    if(max_seg_len < min_seg_len)
    {
        stop("max_seg_len must be greater than min_seg_len")
    }
    # check beta_tilde
    if(!is.null(beta_tilde))
        {
            if(!is_numeric(beta_tilde))
            {
                stop("beta_tilde must be numeric")
            }
            if(length(beta_tilde) != 1)
            {
                stop("beta_tilde must be a single scalar value")
            }
            if(beta_tilde < 0)
            {
                stop("beta_tilde must be >= 0")
            }
        }
    # check beta
    if(!is.null(beta))
    {
        if(!is_numeric(beta))
        {
            stop("beta must be numeric")
        }
        if(!(all(beta >= 0)))
        {
                stop("beta values must be >= 0")
        }
        if(!(all(beta_tilde >= 0)))
        {
                stop("beta_tilde values must be >= 0")
        }
    }
    # wrap the callable transform object in a function to store in S4 class
    transform_method<-function(arg){return(transform(arg))}
    # call appropriate helper function
    tryCatch(
    {
        if(univariate)
        {
            res<-capa.uv_call(x,beta,beta_tilde,type,min_seg_len,max_seg_len)
            res@data<-x_untransformed
            res@transform<-transform_method
            return(res)
        }
        else
        {
            res<-capa.mv_call(x,beta,beta_tilde,type,min_seg_len,max_seg_len,max_lag)
            res@data<-x_untransformed
            res@transform<-transform_method
            return(res)
        }
    },error = function(e) {print(e$message);stop();})
    
}


capa_line_plot<-function(object,epoch=dim(object@data)[1],subset=1:ncol(object@data),variate_names=FALSE)
{
    # creating null entries for ggplot global variables so as to pass CRAN checks
    x<-value<-ymin<-ymax<-x1<-x2<-y1<-y2<-x1<-x2<-y1<-y2<-NULL
    data_df<-as.data.frame(object@data)
    names<-paste("y",1:ncol(object@data),sep="")
    colnames(data_df)<-names
    data_df<-as.data.frame(data_df[,subset,drop=FALSE])
    n<-nrow(data_df)
    p<-ncol(data_df)
    data_df<-cbind(data.frame("x"=1:n),data_df)
    data_df<-melt(data_df,id="x")
    out<-ggplot(data=data_df)
    out<-out+aes(x=x,y=value)
    out<-out+geom_point(alpha=0.3)
    # highlight the collective anomalies
    c_anoms<-collective_anomalies(object,epoch=epoch)
    c_anoms<-c_anoms[c_anoms$variate %in% subset,]
    if(!any(is.na(c_anoms)) & nrow(c_anoms) > 0)
    {
        c_anoms_data_df<-c_anoms[,1:3]
        names(c_anoms_data_df)<-c(names(c_anoms_data_df)[1:2],"variable")
        c_anoms_data_df$variable<-names[c_anoms_data_df$variable]
        c_anoms_data_df$ymin<--Inf
        c_anoms_data_df$ymax<-Inf
        out<-out+geom_rect(data=c_anoms_data_df,inherit.aes = F,mapping=aes(xmin=start,xmax=end,ymin=ymin,ymax=ymax),fill="blue",alpha=0.2)
        c_anoms_data_df$start<-c_anoms_data_df$start+c_anoms$start.lag
        c_anoms_data_df$end<-c_anoms_data_df$end-c_anoms$end.lag
        out<-out+geom_rect(data=c_anoms_data_df,inherit.aes = F,mapping=aes(xmin=start,xmax=end,ymin=ymin,ymax=ymax),fill="blue",alpha=0.5)
    }
    # out<-out+facet_grid(variable~.,scales="free_y")
    # highlight the point anomalies
    p_anoms<-point_anomalies(object,epoch)
    p_anoms<-p_anoms[p_anoms$variate %in% subset,]
    if(!any(is.na(p_anoms)) & nrow(p_anoms) > 0)
        {
            p_anoms_data_df<-Reduce(rbind,Map(function(a,b) data_df[data_df$variable==names[a] & data_df$x==b,],p_anoms$variate,p_anoms$location))
            out<-out+geom_point(data=p_anoms_data_df,colour="red", size=1.5,alpha=0.3)
        }
    out<-out+facet_grid(factor(variable,levels=(names)) ~ .,scales="free_y")
    # grey out the data after epoch
    if(epoch != nrow(object@data))
        {
	    d<-data.frame(variable=names[subset],x1=epoch,x2=n,y1=-Inf,y2=Inf)	
            out<-out+geom_rect(data=d,inherit.aes=F,mapping=aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2),fill="yellow",alpha=0.2)
        }
    out<-out+theme_bw()
    out<-out+theme(axis.text.y=element_blank())
    out<-out+labs(x="t")
    if(variate_names==FALSE)
        {
            out<-out+theme(strip.text.y=element_blank())
        }
    return(out)
}


capa_tile_plot<-function(object,variate_names=FALSE,epoch=dim(object@data)[1],subset=1:ncol(object@data))
{
    # nulling out variables used in ggplot to get the package past CRAN checks
    x1<-y1<-x2<-y2<-variable<-value<-NULL
    df<-as.data.frame(object@data)
    df<-as.data.frame(df[,rev(subset),drop=FALSE])
    # normalise data
    for(i in 1:ncol(df))
    {
        df[,i]<-(df[,i]-min(df[,i]))/(max(df[,i])-min(df[,i]))
    }
    t<-data.frame("t"=seq(1,nrow(df)))
    molten.data<-melt(cbind(t,df),id="t")
    out<-ggplot(molten.data, aes(t,variable))
    out<-out+geom_tile(aes(fill=value))
    # get any collective anomalies
    c_anoms<-collective_anomalies(object,epoch=epoch)
    c_anoms<-c_anoms[c_anoms$variate %in% subset,]
    c_anoms<-unique(c_anoms[,1:2])
    if(!any(is.na(c_anoms)) & nrow(c_anoms) > 0)
    {
        ymin<-0
        ymax<-ncol(df)
        out<-out+annotate("rect",xmin=c_anoms$start,xmax=c_anoms$end,ymin=ymin,ymax=ymax+1,alpha=0.0,color="red",fill="yellow")
    }
    if(epoch != nrow(object@data))
        {
            d<-data.frame(x1=epoch,x2=nrow(object@data),y1=-Inf,y2=Inf)
            out<-out+geom_rect(data=d,inherit.aes=F,mapping=aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2),fill="yellow",alpha=0.2)
        }

    out<-out+theme_bw()
    out<-out+theme(axis.text.y=element_blank())
    out<-out+theme(axis.ticks.y=element_blank())
    out<-out+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())


    return(out)
}


#' Visualisation of data, collective and point anomalies.
#'
#' @name plot
#'
#' @description Plot methods for S4 objects returned by \code{\link{capa}}, \code{\link{capa.uv}}, \code{\link{capa.mv}}, \code{\link{scapa.uv}}, \code{\link{scapa.mv}}, \code{\link{pass}}, and \code{\link{sampler}}. 
#'
#' The plot can either be a line plot or a tile plot, the type produced depending on the options provided to the \code{plot} function and/or the dimensions of the
#' data associated with the S4 object.
#'
#' @docType methods
#'
#' @param x An instance of an S4 class produced by \code{\link{capa}}, \code{\link{capa.uv}}, \code{\link{capa.mv}}, \code{\link{pass}}, or \code{\link{sampler}}. 
#' @param subset A numeric vector specifying a subset of the variates to be displayed. Default value is all of the variates present in the data.
#' @param variate_names Logical value indicating if variate names should be displayed on the plot. This is useful when a large number of variates are being displayed
#' as it makes the visualisation easier to interpret. Default value is FALSE.
#' @param tile_plot Logical value. If TRUE then a tile plot of the data is produced. The data displayed in the tile plot is normalised to values in [0,1] for each variate.
#' This type of plot is useful when the data contains are large number of variates. The default value is TRUE if the number of variates is greater than 20.
#' @param epoch Positive integer. CAPA methods are sequential and as such, can generate results up to, and including, any epoch within the data series. This can be controlled by the value
#' of \code{epoch} and is useful for examining how the inferred anomalies are modified as the data series grows. The default value for \code{epoch} is the length of the data series.
#' 
#' @return A ggplot object.
#'
#' @rdname plot-methods
#'
#' @aliases plot,capa.class-method
#' 
#' @seealso \code{\link{capa}},\code{\link{capa.uv}},\code{\link{capa.mv}},\code{\link{pass}},\code{\link{sampler}}.
#'
#' @export 
setMethod("plot",signature=list("capa.class"),function(x,subset,variate_names=FALSE,tile_plot,epoch=nrow(x@data))
{
    if(missing(subset))
    {
        subset<-1:ncol(x@data)
    }
    if(epoch < 0)
    {
        stop("epoch should be a positive integer")
    }
    if(epoch > nrow(x@data))
    {
        stop("epoch cannot be greater than the number of observations in the data")
    }
    if(missing(tile_plot))
    {
        tile_plot<-NULL
    }
    subset<-sort(unique(subset))
    if(!is.logical(tile_plot))
    {
        if(is.null(tile_plot))
        {
            tile_plot<-FALSE
            if(length(subset) > 20)
            {
                tile_plot<-TRUE
            }
        }
        else
        {
            stop("tile_plot must be of type logical or NULL")
        }
    }
    if(!is.logical(variate_names))
    {
       stop("variable_names must be of type logical or NULL")
    }
    if(tile_plot)
    {
        return(capa_tile_plot(x,variate_names=FALSE,epoch=epoch,subset=subset))
    }
    else
    {
        return(capa_line_plot(x,variate_names=variate_names,epoch=epoch,subset=subset))
    }
})







