#' Running mean of a vector
#' 
#' Function to calculate the running mean of a numeric vector
#' 
#' 
#' @param vec numeric vector
#' @param runn_mean number of vector elements to use for calculating the
#' running mean
#' @param na.rm ignore NA values when calculating means. Defaults to FALSE.
#' @param exclude_central_value exclude central value in calculating means.
#' Defaults to FALSE.
#' @param FUN function to be applied. For a running mean, this is usually mean (the
#' default), but other functions can also be specified here (the na.rm parameter
#' won't work then, and the function has to be dependent on one numeric variable only.
#' @return numeric vector containing the running mean
#' @author Eike Luedeling
#' @keywords "running mean"
#' @examples
#' 
#' 
#' plot(runn_mean(rnorm(1000),150))
#' 
#' @export runn_mean
runn_mean<-function(vec,runn_mean,na.rm=FALSE,exclude_central_value=FALSE,FUN=mean)
{
  
  runn_mean<-min(c(runn_mean,floor(length(vec)/2)))
  
  if(identical(FUN,mean)) if(na.rm==TRUE) FUN<-function(x) mean(x,na.rm=TRUE) else FUN<-function(x) mean(x,na.rm=FALSE)
  
  ww <- vec
  rr <- vec
  for (dd in 1:length(ww)) {
    if (dd < ceiling(runn_mean/2)) {
      if(!exclude_central_value) rr[dd] <- sapply(list(ww[1:(dd + floor(runn_mean/2))]),FUN)
      if(exclude_central_value) rr[dd] <- sapply(list(ww[(1:(dd + floor(runn_mean/2)))[which(!(1:(dd + floor(runn_mean/2)))==dd)]]),FUN)
    }
    if ((dd >= ceiling(runn_mean/2)) & (dd <= length(ww) - 
                                        ceiling(runn_mean/2))) {
      if(!exclude_central_value) rr[dd] <- sapply(list(ww[(dd - floor(runn_mean/2)):(dd + 
                                                                                       floor(runn_mean/2))]),FUN)
      if(exclude_central_value) rr[dd] <- sapply(list(ww[((dd - floor(runn_mean/2)):(dd + 
                                                                                       floor(runn_mean/2)))[
                                                                                         which(!((dd - floor(runn_mean/2)):(dd + 
                                                                                                                              floor(runn_mean/2)))==dd)]]),FUN)
    }
    if (dd > (length(ww) - ceiling(runn_mean/2))) {
      if(!exclude_central_value) rr[dd] <- sapply(list(ww[(dd - floor(runn_mean/2)):length(ww)]),FUN)
      if(exclude_central_value) rr[dd] <- sapply(list(ww[((dd - floor(runn_mean/2)):length(ww))[
        which(!((dd - floor(runn_mean/2)):length(ww))==dd)]]),FUN)
    }
  }
  return(rr)}


#' Prediction based on a running mean
#' 
#' Function to predict values based on a running mean (or another function) of a numeric vector.
#' 
#' 
#' @param indep numeric vector of independent variables, should be sequential
#' @param dep numeric vector of dependent variables
#' @param pred numeric vector of values to be predicted
#' @param runn_mean number of vector elements to use for calculating the
#' running mean
#' @param na.rm ignore NA values when calculating means. Defaults to FALSE.
#' @param exclude_central_value exclude central value in calculating means.
#' Defaults to FALSE.
#' @param FUN function to be applied. For a running mean, this is usually mean (the
#' default), but other functions can also be specified here (the na.rm parameter
#' won't work then, and the function has to be dependent on one numeric variable only.
#' 
#' @details The running mean calculation that underlies the prediction is based purely on the sequence of
#' observed values, without accounting for any variation in intervals of the independent data. This
#' means that the function performs best with regularly spaced independent variables. Note that the
#' function will return NA when asked to predict values that are outside the range of independent values
#' provided as input. The prediction results are computed by linearly interpolating between the running
#' mean values determined for the nearest neighbors of the value that is to be predicted.
#' 
#' @return list of two elements, with $x containing the values to be predicted and $predicted
#' the predicted values
#' @author Eike Luedeling
#' @keywords "running mean" prediction
#' @examples
#' 
#' indep<-(1:100)
#' dep<-sin(indep/20)+rnorm(100)/5
#' pred<-c(12,13,51,70,90)
#' 
#' predicted<-runn_mean_pred(indep,dep,pred,runn_mean = 25)
#' 
#' plot(dep~indep)
#' points(predicted$predicted~predicted$x,col="red",pch=15)
#' 
#' @export runn_mean_pred
runn_mean_pred<-function(indep,dep,pred,runn_mean=11,na.rm=FALSE,exclude_central_value=FALSE,FUN=mean)
{
  
  #if(!(min(indep,na.rm=TRUE)<min(pred,na.rm=TRUE)&max(pred,na.rm=TRUE)<max(indep,na.rm=TRUE))) stop("Error: Value to be predicted is outside the range defined by the data")
  
  runn_mean<-min(c(runn_mean,floor(length(indep)/2)))
  
  runny<-runn_mean(dep,runn_mean,na.rm=na.rm,exclude_central_value=exclude_central_value,FUN=FUN)
  
  pred_fun<-function(pr)
  {
    if(is.na(pr)) return(NA)
    if(!(min(indep,na.rm=TRUE)<=pr&pr<=max(indep,na.rm=TRUE))) return(NA)
    
    start_y<-runny[max(which(indep<=pr))]
    end_y<-runny[min(which(indep>=pr))]
    
    start_x<-indep[max(which(indep<=pr))]
    end_x<-indep[min(which(indep>=pr))]
    
    if(end_x==start_x) out<-runny[which(indep==pr)] else
      out<-start_y+(end_y-start_y)*((pr-start_x)/(end_x-start_x))
    return(out)
  }
  
  return(list(x=pred,predicted=sapply(pred,pred_fun)))
}

