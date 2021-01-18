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
#' @include global.R kmmNames.R IClusterModel.R
NULL

#-----------------------------------------------------------------------
#' Create an instance of the [\code{\linkS4class{KmmMixedDataModel}}] class
#'
#' This function computes the optimal mixture model for mixed data using kernel
#' mixture models according to the \code{criterion} among the number of clusters
#' given in \code{nbCluster} using the strategy specified in [\code{strategy}].
#'
#' @param ldata [\code{list}] containing the data sets (matrices and/or data.frames).
#' @param lmodels a [\code{list}] of same length than data. It contains the model
#' names, kernel names and kernel parameter names to use in order to fit each
#' data set.
#' @param nbCluster [\code{\link{vector}}] with the number of clusters to test.
#' @param strategy a [\code{\linkS4class{ClusterStrategy}}] object containing
#' the strategy to run. Default is clusterStrategy().
#' @param criterion character defining the criterion to select the best model.
#' The best model is the one with the lowest criterion value.
#' Possible values: "BIC", "AIC", "ICL", "ML". Default is "ICL".
#' @param nbCore integer defining the number of processors to use (default is 1, 0 for all).
#'
#' @details For each data set in data, we need to specify a list of parameters
#' 
#' @examples
#'
#' ## An example with the bullsEye data set
#' data(bullsEye)
#' data(bullsEye.cat)
#' ## with default values
#' ldata     <- list(bullsEye, bullsEye.cat)
#' modelcont <- list(modelName="kmm_pk_s", dim = 10, kernelName="Gaussian")
#' modelcat  <- list(modelName="kmm_pk_s", dim = 20, kernelName="Hamming", kernelParameters = c(0.6))
#' lmodels   <- list( modelcont, modelcat)
#' 
#' model <- kmmMixedData(ldata, lmodels, nbCluster=2:5, strategy = clusterFastStrategy())
#'
#' ## get summary
#' summary(model)
#'
#'
#' \dontrun{
#' ## use graphics functions
#' plot(model)
#' }
#'
#' @return An instance of the [\code{\linkS4class{KmmMixedDataModel}}] class.
#' @author Serge Iovleff
#'
kmmMixedData <- function( ldata, lmodels, nbCluster=2
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
  { stop("criterion is not valid. See ?KmmMixedDataModel for the list of valid criterion")}
  
  # check strategy
  if(class(strategy)[1] != "ClusterStrategy")
  {stop("strategy is not a Cluster Stategy class (must be an instance of the class ClusterStrategy).")}
  validObject(strategy);
  
  # check ldata and lmodels
  if (!is.list(ldata))   { stop("data must be a list")}
  if (!is.list(lmodels)) { stop("lmodels must be a list")}
  if (length(ldata) != length(lmodels)) { stop("data and lmodels must be of equal lengths");}
  
  # create list of models
  lcomponent <- vector("list", length(ldata));
  for (i in 1:length(ldata))
  {
    # check data
    ldata[[i]] <- as.matrix(ldata[[i]])
    # check parameter
    param <- lmodels[[i]];
    if (is.list(param))
    {
      modelName         <- param$modelName
      dim               <- param$dim
      kernelName        <- param$kernelName
      kernelParameters  <- param$kernelParameters
      kernelComputation <- param$kernelComputation
      # check and set default values
      if (is.null(kernelName))        { kernelName <- "Gaussian"}
      if (is.null(dim))               { dim <- 10}
      if (is.null(kernelParameters))  { kernelParameters <- c(1)}
      if (is.null(kernelComputation)) { kernelComputation <- TRUE}
    }
    else
    { 
      modelName  <- param
      dim        <- 10
      kernelName <- "Gaussian"
      kernelParameters <- c(1)
      kernelComputation <- TRUE
    }
    if (!kmmValidModelNames(modelName))
    { stop("modelName is not valid. See ?kmmNames for the list of valid model names")}
    if (!kmmValidKernelNames(kernelName))
    { stop("kernelName is not valid. See ?kmm for the list of valid model names")}
    
    # create component
    lcomponent[[i]] = new( "KmmComponent", ldata[[i]], dim
                         , nbClusterMin
                         , modelName
                         , kernelName, kernelParameters, kernelComputation)
  } # for i
  
  # Create model
  model = new("KmmMixedDataModel", lcomponent)
  model@strategy      = strategy;
  model@criterionName = criterion
  
  # start estimation of the models
  resFlag  <- FALSE;
  if (length(nbCluster) >0)
  {
    resFlag = .Call("kmmMixedData", model, nbCluster, nbCore, PACKAGE="MixAll");
  }
  # set names
  if (resFlag != TRUE) {cat("WARNING: An error occurs during the clustering process");}
  model
}

#-----------------------------------------------------------------------
#' Definition of the [\code{\linkS4class{KmmMixedDataModel}}] class
#'
#' This class defines a mixed data kernel mixture Model (KMM).
#'
#' This class inherits from the [\code{\linkS4class{IClusterModel}}] class.
#' A model for mixed data is a mixture model of the form:
#' \deqn{
#' f({{x}}_i=({{x}}_{1i}, {{x}}_{2i},\ldots {{x}}_{Li})|\theta)
#' = \sum_{k=1}^K p_k \prod_{l=1}^L h({{x}}_{li}).
#' }
#' The density functions (or probability distribution functions)
#' \deqn{h(.)} can be any implemented kmm model on a RKHS space.
#'
#' @slot lcomponent  a list of [\code{\linkS4class{KmmComponent}}]
#' @seealso [\code{\linkS4class{IClusterModel}}] class
#'
#' @examples
#' getSlots("KmmMixedDataModel")
#'
#' @author Serge Iovleff
#'
#' @name KmmMixedDataModel
#' @rdname KmmMixedDataModel-class
#' @aliases KmmMixedDataModel-class
#'
setClass(
  Class = "KmmMixedDataModel",
  representation( lcomponent = "list"),
  contains=c("IClusterModel"),
  validity=function(object)
  {
    nbData = length(object@lcomponent[1])
    if (nbData == 0) { stop("At least on data set must be given.")}
    for (l in 1:nbData)
    {
      if (nrow(object@lcomponent[[1]]@data) != object@nbSample)
      { stop("All data sets must have the same number of individuals (number of rows)")}
    }
    return(TRUE)
  }
)

#' Initialize an instance of a MixAll S4 class.
#'
#' Initialization method of the [\code{\linkS4class{KmmMixedDataModel}}] class.
#' Used internally in the 'MixAll' package.
#'
#' @rdname initialize-methods
#' @keywords internal
setMethod(
  f="initialize",
  signature=c("KmmMixedDataModel"),
  definition=function(.Object, lcomponent, nbCluster=2)
  {
    # for data
    if(missing(lcomponent)) {stop("lcomponent is mandatory in KmmMixedDataModel.")}
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
#' @aliases print print,KmmMixedDataModel-method
setMethod(
  f="print",
  signature=c("KmmMixedDataModel"),
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
#' @aliases show-KmmMixedDataModel,KmmMixedDataModel,KmmMixedDataModel-method
setMethod(
  f="show",
  signature=c("KmmMixedDataModel"),
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
#' @aliases summary summary,KmmMixedDataModel-method
setMethod(
  f="summary",
  signature=c("KmmMixedDataModel"),
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

#' Plotting of a class [\code{\linkS4class{KmmMixedDataModel}}]
#'
#' Plotting data from a [\code{\linkS4class{KmmMixedDataModel}}] object
#' using the estimated parameters and partition.
#'
#' @param x an object of class [\code{\linkS4class{KmmMixedDataModel}}]
#' @param y a vector listing the data sets you want to disply
#' @param ... further arguments passed to or from other methods
#'
#' @aliases plot-KmmMixedDataModel
#' @docType methods
#' @rdname plot-KmmMixedDataModel-method
#'
#'
#' @seealso \code{\link{plot}}
#' @examples
#' \dontrun{
#' data(bullsEye)
#' data(bullsEye.cat)
#' ## with default values
#' ldata  = list(bullsEye, bullsEye.cat)
#' modelcont <- list(modelName="kmm_pk_s", dim = 10, kernelName="Gaussian")
#' modelcat  <- list(modelName="kmm_pk_s", dim = 20, kernelName="Hamming", kernelParameters = c(0.6))
#' lmodels = list( modelcont, modelcat)
#' 
#' model <- kmmMixedData(ldata, lmodels, nbCluster=2:5, strategy = clusterFastStrategy())
#' # plot only the first continuous data set
#' plot(model, y=c(1))
#'   }
#'
setMethod(
  f="plot",
  signature=c("KmmMixedDataModel"),
  function(x, y, ...)
  {
    # total number of cluster in the data set
    nbCluster   <- ncol(x@tik);
    nbComponent <- length(x@lcomponent)
    #
    # check y, no y => display all dimensions
    if (missing(y)) { y=1:nbComponent; }
    else
    {
      if (max(y)>nbComponent)
        stop("y should not be greater than the number of data set")
    }
    #
    par(mfrow=c(1, length(y)), pty="s")
    for (i in y)
    {
      plot(x@lcomponent[[i]], x@zi,...)
    }
  }
)

