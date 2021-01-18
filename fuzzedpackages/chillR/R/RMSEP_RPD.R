#' Root Mean Square Error of Prediction (RMSEP)
#' 
#' This function computes the Root Mean Square Error of Prediction (RMSEP),
#' a commonly used measure for the predictive capacity of a model. It
#' compares values predicted with a model with observed values.
#' 
#' @param predicted a numeric vector containing predicted values.
#' @param observed a numeric vector of the same length as ```predicted```
#' containing observed values.
#' @param na.rm Boolean parameter indicating whether NA values should be removed before the analysis
#' @return numeric value of the RMSEP.
#' @author Eike Luedeling
#' @keywords model validation
#' @examples
#' 
#' predicted<-c(1,2,3,4,5,6,7,8,9,10)
#' observed<-c(1.5,1.8,3.3,3.9,4.4,6,7.5,9,11,10)
#' 
#' RMSEP(predicted,observed)
#' 
#' @export RMSEP
RMSEP<-function(predicted,observed,na.rm=FALSE)
{if(!na.rm)
  if(!(length(which(is.na(predicted)))+length(which(is.na(observed))))==0)
    stop("Datasets include NA values. This may indicate a serious prediction problem. To override this error, set na.rm=TRUE.")
  sqrt(sum((observed-predicted)^2,na.rm=TRUE)/
         length(which(!is.na(observed-predicted))))}


#' Residual Prediction Deviation (RPD)
#' 
#' This function computes the Residual Prediction Deviation (RPD), which
#' is defined as the standard deviation of observed values divided by
#' the Root Mean Square Error or Prediction (RMSEP). The RDP takes both
#' the prediction error and the variation of observed values into 
#' account, providing a metric of model validity that is more objective
#' than the RMSEP and more easily comparable across model validation
#' studies. The greater the RPD, the better the model's predictive
#' capacity.
#' 
#' Interpretation of the RPD is somewhat arbitrary, with different
#' thresholds for a good model used in the literature. Many studies
#' call a model *excellent*, when the RPD is above 2 (but other
#' classification use thresholds as high as 8 for this).
#' 
#' @param predicted a numeric vector containing predicted values.
#' @param observed a numeric vector of the same length as ```predicted```
#' containing observed values.
#' @param na.rm Boolean parameter indicating whether NA values should be removed before the analysis
#' @return numeric value of the RDP.
#' @author Eike Luedeling
#' @keywords model validation
#' @references Williams PC and Sobering DC (1993) Comparison of
#' commercial near infrared transmittance and reflectance instruments
#' for analysis of whole grains and seeds. J. Near Infrared Spectrosc.
#' 1, 25-32 (I didn't have access to this paper, but have noticed that
#' it is often provided as the key reference for the RPD).
#' 
#' @examples
#' 
#' predicted<-c(1,2,3,4,5,6,7,8,9,10)
#' observed<-c(1.5,1.8,3.3,3.9,4.4,6,7.5,9,11,10)
#' 
#' RPD(predicted,observed)
#' 
#' @export RPD
RPD<-function(predicted,observed,na.rm=FALSE)
{if(!na.rm)
  if(!(length(which(is.na(predicted)))+length(which(is.na(observed))))==0)
    stop("Datasets include NA values. This may indicate a serious prediction problem. To override this error, set na.rm=TRUE.")
  sd(observed,na.rm=TRUE)/
    sqrt(sum((observed-predicted)^2,na.rm=TRUE)/
           length(which(!is.na(observed-predicted))))}

#' Ratio of Performance to InterQuartile distance (RPIQ)
#' 
#' This function computes the Ratio of Performance to InterQuartile
#' distance (RPIQ), which
#' is defined as interquartile range of the observed values divided by
#' the Root Mean Square Error or Prediction (RMSEP). The RPIQ takes both
#' the prediction error and the variation of observed values into 
#' account, providing a metric of model validity that is more objective
#' than the RMSEP and more easily comparable across model validation
#' studies. The greater the RPIQ, the better the model's predictive
#' capacity. In contrast to the Residual Prediction Deviation (RPD),
#' the RPIQ makes no assumptions about the distribution of the
#' observed values (since the RDP includes a standard deviation, it
#' assumed normal distribution of the observed values).
#' 
#' Interpretation of the RPIQ differs in the literature, with different
#' thresholds used for judging model quality.
#' 
#' @param predicted a numeric vector containing predicted values.
#' @param observed a numeric vector of the same length as ```predicted```
#' containing observed values.
#' @param na.rm Boolean parameter indicating whether NA values should be removed before the analysis
#' @return numeric value of the RPIQ
#' @author Eike Luedeling
#' @keywords model validation
#' @references Bellon-Maurel V, Fernandez-Ahumada E, Palagos B,
#' Roger J-M, McBratney A, 2010. Critical review of chemometric
#' indicators commonly used for assessing the quality of the prediction
#' of soil attributes by NIR spectroscopy, In TrAC Trends in Analytical
#' Chemistry 29(9), 1073-1081.
#' @importFrom stats quantile
#' 
#' @examples
#' 
#' predicted<-c(1,2,3,4,5,6,7,8,9,10)
#' observed<-c(1.5,1.8,3.3,3.9,4.4,6,7.5,9,11,10)
#' 
#' RPD(predicted,observed)
#' 
#' @export RPIQ
RPIQ<-function(predicted,observed,na.rm=FALSE)
{if(!na.rm)
  if(!(length(which(is.na(predicted)))+length(which(is.na(observed))))==0)
    stop("Datasets include NA values. This may indicate a serious prediction problem. To override this error, set na.rm=TRUE.")
  as.numeric(quantile(observed,na.rm=na.rm)[4]-quantile(observed,na.rm=na.rm)[2])/
    sqrt(sum((observed-predicted)^2,na.rm=na.rm)/
           length(which(!is.na(observed-predicted))))}
