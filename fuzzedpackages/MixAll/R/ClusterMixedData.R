#-----------------------------------------------------------------------
#     Copyright (C) 2012-2016  Serge Iovleff, University Lille 1, Inria
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as
#    published by the Free Software Foundation; either version 2 of the
#    License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public
#    License along with this program; if not, write to the
#    Free Software Foundation, Inc.,
#    59 Temple Place,
#    Suite 330,
#    Boston, MA 02111-1307
#    USA
#
#    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
#
#-----------------------------------------------------------------------
#' @include global.R ClusterModelNames.R IClusterModel.R
NULL

#-----------------------------------------------------------------------
#' Create an instance of the [\code{\linkS4class{ClusterMixedDataModel}}] class
#'
#' This function computes the optimal mixture model for mixed data according
#' to the \code{criterion} among the number of clusters given in
#' \code{nbCluster} using the strategy specified in [\code{strategy}].
#'
#' @param data [\code{list}] containing the data sets (matrices and/or data.frames).
#' If data sets contain NA values, these missing values will be estimated during
#' the estimation process.
#' @param models a [\code{vector}] of character or a [\code{list}] of
#' same length than data. It contains the model names to use in order to fit
#' each data set.
#' @param nbCluster [\code{\link{vector}}] with the number of clusters to test.
#' @param strategy a [\code{\linkS4class{ClusterStrategy}}] object containing
#' the strategy to run. Default is clusterStrategy().
#' @param criterion character defining the criterion to select the best model.
#' The best model is the one with the lowest criterion value.
#' Possible values: "BIC", "AIC", "ICL", "ML". Default is "ICL".
#' @param nbCore integer defining the number of processors to use (default is 1, 0 for all).
#'
#' @examples
#' ## A quantitative example with the heart disease data set
#' data(HeartDisease.cat)
#' data(HeartDisease.cont)
#' ## with default values
#' ldata = list(HeartDisease.cat, HeartDisease.cont);
#' models = c("categorical_pk_pjk","gaussian_pk_sjk")
#' model <- clusterMixedData(ldata, models, nbCluster=2:5, strategy = clusterFastStrategy())
#'
#' ## get summary
#' summary(model)
#'
#' ## get estimated missing values
#' missingValues(model)
#'
#' \dontrun{
#' ## print model
#' print(model)
#' ## use graphics functions
#' plot(model)
#' }
#'
#' @return An instance of the [\code{\linkS4class{ClusterMixedDataModel}}] class.
#' @author Serge Iovleff
#'
clusterMixedData <- function( data, models, nbCluster=2
                            , strategy=clusterStrategy()
                            , criterion="ICL"
                            , nbCore = 1)
{
  # check nbCluster
  nbClusterModel = length(nbCluster);
  nbClusterMin   = min(nbCluster);
  nbClusterMax   = max(nbCluster);
  if (nbClusterMin < 1) { stop("The number of clusters must be greater or equal to 1")}
  # check criterion
  if(sum(criterion %in% c("BIC","AIC", "ICL", "ML")) != 1)
  { stop("criterion is not valid. See ?clusterMixedData for the list of valid criterion")}
  # check strategy
  if(class(strategy)[1] != "ClusterStrategy")
  {stop("strategy is not a Cluster Stategy class (must be an instance of the class ClusterStrategy).")}
  validObject(strategy);
  # check data and models
  if (!is.list(data))     { stop("data must be a list");}
  if (!is.vector(models)) { stop("models must be a vector of character");}
  if (length(data) != length(models)) { stop("data and models must be of equal lengths");}
  
  # create list of component
  lcomponent <- vector("list", length(data));
  for (i in 1:length(data))
  {
    param <- models[i];
    modelName <- models[i];
    # check if it is a Categorical model
    if( clusterValidCategoricalNames(modelName) )
    { lcomponent[[i]] <- new("ClusterCategoricalComponent", data[[i]], nbClusterMin, modelName)}
    else
    {  # check if it is a Gamma model
      if( clusterValidGammaNames(modelName) )
      { lcomponent[[i]] <- new("ClusterGammaComponent", data[[i]], nbClusterMin, modelName)}
      else
      { # check if it is a diagonal Gaussian model
        if( clusterValidDiagGaussianNames(modelName) )
        { lcomponent[[i]] <- new("ClusterDiagGaussianComponent", data[[i]], nbClusterMin, modelName);}
        else
        { # check if it is a Poisson model
          if( clusterValidPoissonNames(modelName) )
          { lcomponent[[i]] <- new("ClusterPoissonComponent", data[[i]], nbClusterMin, modelName);}
          else
          {
            stop("in clusterMixedData: invalid model name");
          }
        }
      } # else gamma
    } # else categorical
  } # for i
  # Create model
  model = new("ClusterMixedDataModel", lcomponent)
  model@strategy = strategy;
  model@criterionName = criterion
  
  # start estimation of the models
  resFlag  <- FALSE;
  if (length(nbCluster) >0)
  {
    resFlag = .Call("clusterMixedData", model, nbCluster, nbCore, PACKAGE="MixAll");
  }
  # set names
  if (resFlag != TRUE) {cat("WARNING: An error occurs during the clustering process");}
  for (i in 1:length(data))
  {
    if(clusterValidCategoricalNames(models[i]))
    { dim(model@lcomponent[[i]]@plkj) <- c(model@lcomponent[[i]]@nbModalities, model@nbCluster, ncol(model@lcomponent[[i]]@data))}
  }
  model
}

#-----------------------------------------------------------------------
#' Definition of the [\code{\linkS4class{ClusterMixedDataModel}}] class
#'
#' This class defines a mixed data mixture Model.
#'
#' This class inherits from the [\code{\linkS4class{IClusterModel}}] class.
#' A model for mixed data is a mixture model of the form:
#' \deqn{
#' f({{x}}_i=({{x}}_{1i}, {{x}}_{2i},\ldots {{x}}_{Li})|\theta)
#' = \sum_{k=1}^K p_k \prod_{l=1}^L h({{x}}_{li}| \lambda_{lk},\alpha_l).
#' }
#' The density functions (or probability distribution functions)
#' \deqn{h(.|\lambda_{lk},\alpha_l)}
#' can be any implemented model (Gaussian, Poisson,...).
#'
#' @slot lcomponent  a list of [\code{\linkS4class{IClusterComponent}}]
#' @seealso [\code{\linkS4class{IClusterModel}}] class
#'
#' @examples
#' getSlots("ClusterMixedDataModel")
#'
#' @author Serge Iovleff
#'
#' @name ClusterMixedDataModel
#' @rdname ClusterMixedDataModel-class
#' @aliases ClusterMixedDataModel-class
#'
setClass(
  Class = "ClusterMixedDataModel",
  representation( lcomponent = "list"),
  contains=c("IClusterModel"),
  validity=function(object)
  {
    nbData = length(object@lcomponent)
    if (nbData == 0) {stop("At least on data set must be given.");}
    for (l in 1:nbData)
    {
      if (nrow(object@lcomponent[[1]]@data) != object@nbSample)
      {stop("All data sets must have the same number of individuals (number of rows).");}
    }
    return(TRUE)
  }
)

#' Initialize an instance of a MixAll S4 class.
#'
#' Initialization method of the [\code{\linkS4class{ClusterMixedDataModel}}] class.
#' Used internally in the 'MixAll' package.
#'
#' @rdname initialize-methods
#' @keywords internal
setMethod(
  f="initialize",
  signature=c("ClusterMixedDataModel"),
  definition=function(.Object, lcomponent, nbCluster=2)
  {
    # for data
    if(missing(lcomponent)) {stop("lcomponent is mandatory in ClusterMixedDataModel.")}
    nbData = length(lcomponent)
    if (nbData == 0) {stop("At least on data set must be given.")}
    .Object@lcomponent <- lcomponent;
    # take first element of the list, this will give us the dimensions
    nbSample = nrow(.Object@lcomponent[[1]]@data);
    .Object <- callNextMethod(.Object, nbSample, nbCluster)
    # validate
    validObject(.Object)
    return(.Object)
  }
)

#' @rdname print-methods
#' @aliases print print,ClusterMixedDataModel-method
setMethod(
  f="print",
  signature=c("ClusterMixedDataModel"),
  function(x,...){
    cat("****************************************\n")
    callNextMethod()
    nbData <- length(x@lcomponent)
    if(nbData>0)
    {
      for (l in 1:nbData)
      {
        cat("* model name = ", x@lcomponent[[l]]@modelName, "\n")
        print(format(x@lcomponent[[l]]@data),quote=FALSE);
      }
    }
    cat("****************************************\n")
    if(nbData>0)
    {
      for (l in 1:nbData)
      {
        cat("* model name = ", x@lcomponent[[l]]@modelName, "\n");
        for(k in 1:length(x@pk))
        {
          cat("*** Cluster: ",k,"\n")
          cat("* Proportion = ", format(x@pk[k]), "\n")
          print(x@lcomponent[[l]],k);
        }
      }
      cat("****************************************\n")
    }
  }
)

#' @rdname show-methods
#' @aliases show-ClusterMixedDataModel,ClusterMixedDataModel,ClusterMixedDataModel-method
setMethod(
  f="show",
  signature=c("ClusterMixedDataModel"),
  function(object)
  {
    cat("****************************************\n")
    callNextMethod()
    nbData <- length(object@lcomponent)
    if(nbData>0)
    {
      for (l in 1:nbData)
      {
        cat("* model name = ", object@lcomponent[[l]]@modelName, "\n")
        nrowShow <- min(10,nrow(object@lcomponent[[l]]@data))
        ncolShow <- min(10,ncol(object@lcomponent[[l]]@data))
        cat("* data (limited to 10 samples and 10 variables) =\n")
        print(format(object@lcomponent[[l]]@data[1:nrowShow,1:ncolShow]),quote=FALSE)
      }
    }
    cat("* ... ...\n")
    
    cat("****************************************\n")
    if(nbData>0)
    {
      for (l in 1:nbData)
      {
        cat("* model name = ", object@lcomponent[[l]]@modelName, "\n");
        for(k in 1:length(object@pk))
        {
          cat("*** Cluster: ",k,"\n")
          cat("* Proportion = ", format(object@pk[k]),"\n")
          print(object@lcomponent[[l]], k);
        }
      }
      cat("****************************************\n")
    }
  }
)

#' @rdname summary-methods
#' @aliases summary summary,ClusterMixedDataModel-method
setMethod(
  f="summary",
  signature=c("ClusterMixedDataModel"),
  function(object, ...)
  {
    cat("**************************************************************\n")
    nbData <- length(object@lcomponent)
    if(nbData>0)
    {
      for (l in 1:nbData)
      {
        cat("* model name = ", object@lcomponent[[l]]@modelName, "\n")
      }
    }
    callNextMethod()
    cat("**************************************************************\n")
  }
)

#' Plotting of a class [\code{\linkS4class{ClusterMixedDataModel}}]
#'
#' Plotting data from a [\code{\linkS4class{ClusterMixedDataModel}}] object
#' using the estimated parameters and partition.
#'
#' @param x an object of class [\code{\linkS4class{ClusterMixedDataModel}}]
#' @param y a number between 1 and K-1.
#' @param ... further arguments passed to or from other methods
#'
#' @aliases plot-ClusterMixedDataModel
#' @docType methods
#' @rdname plot-ClusterMixedDataModel-method
#'
#'
#' @seealso \code{\link{plot}}
#' @examples
#' \dontrun{
#'   ## the car data set
#'   data(car)
#'   model <- clusterMixedData(car, 3, strategy = clusterFastStrategy())
#'   plot(model)
#'   }
#'
setMethod(
  f="plot",
  signature=c("ClusterMixedDataModel"),
  function(x, y, ...)
  {
    # total number of cluster in the data set
    nbCluster = ncol(x@tik);
    # check y, no y => display all dimensions
    if (missing(y)) { y=1:(nbCluster-1); }
    else
    { if (round(y)!=y) {stop("y must be an integer.")}
      if (y>nbCluster-1)
        stop("y should not be greater than K-1")
      y <- 1:y
    }
    # get representation
    Y=.visut(x@tik, nbCluster);
    if (nbCluster == 2) { ndim = 1;}
    else { ndim = ncol(Y);}
    # Compute gaussian statistics
    mean  <- matrix(0, nrow = x@nbCluster, ncol =ndim)
    sigma <- matrix(1, nrow = x@nbCluster, ncol =ndim)
    for (k in 1:nbCluster)
    {
      wcov = cov.wt(as.matrix(Y), x@tik[,k], method = "ML");
      mean[k,]  = wcov$center;
      sigma[k,] = sqrt(diag(wcov$cov))
    }
    # create gaussian model
    gauss<-new("ClusterDiagGaussian", Y, nbCluster = x@nbCluster)
    gauss@component@mean = mean
    gauss@component@sigma= sigma
    gauss@pk   = x@pk
    gauss@tik  = x@tik
    gauss@lnFi = x@lnFi
    gauss@zi   = x@zi
    #gauss@component@missing     = x@component@missing
    gauss@lnLikelihood = x@lnLikelihood
    gauss@criterion    = x@criterion
    gauss@nbFreeParameter = x@nbFreeParameter
    gauss@strategy        = x@strategy
    .clusterPlot(gauss, y, .dGauss,...);
  }
)

# get logistic representation
.visut <- function(t, gp)
{ m <- min(t[,gp])
  if (m==0) { t[,gp] = t[,gp] + 1e-30 }
  return(scale(log(sweep(t,1,t[,gp],FUN="/")+ 1e-30), center=TRUE, scale=FALSE)[,-gp])
}
