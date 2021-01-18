#-----------------------------------------------------------------------
#     Copyright (C) 2012-2018  Serge Iovleff, University Lille 1, Inria
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
#' @include IClusterPredict.R
NULL

#-----------------------------------------------------------------------
#' Create an instance of [\code{\linkS4class{ClusterPredict}}] class
#'
#' This function predicts the best cluster each sample in data belongs to. 
#'
#' @param data dataframe or matrix containing the data. Rows correspond to observations
#' and columns correspond to variables. If the data set contains NA values, they
#' will be estimated during the predicting process.
#' @param model (estimated) clustering model to use, i.e. an instance of
#' \code{\linkS4class{ClusterCategorical}}, \code{\linkS4class{ClusterDiagGaussian}},..
#' produced by \code{\link{clusterCategorical}}, \code{\link{clusterDiagGaussian}},...
#' \code{\link{learnCategorical}}, \code{\link{learnDiagGaussian}}, etc.
#' functions.
#' @param algo an instance of \code{\linkS4class{ClusterAlgoPredict}} S4 class. Will not
#' be used if there is no missing values.
#' @param nbCore integer defining the number of processors to use (default is 1, 0 for all).
#'
#' @examples
#' ## A quantitative example with the famous iris data set
#' data(iris)
#' ## get quantitatives 
#' x = as.matrix(iris[1:4])
#' ## sample train and test data sets
#' indexes <- sample(1:nrow(x), nrow(x)/2)
#' train <- x[ indexes,]
#' test  <- x[-indexes,]
#' ## estimate model (using fast strategy, results may be misleading)
#' model1 <- clusterDiagGaussian( data =train, nbCluster=2:3
#'                              , models=c( "gaussian_p_sjk")
#'                              )
#' ## get summary
#' summary(model1)
#' ## compute prediction and compare
#' model2 <- clusterPredict(test, model1)
#' show(model2)
#' as.integer(iris$Species[-indexes])
#'
#' @return An instance of [\code{\linkS4class{ClusterPredict}}] with predicted
#' values
#' @author Serge Iovleff
#'
#'
clusterPredict <- function( data, model, algo = clusterAlgoPredict(), nbCore = 1)
{
  # for model
  if(missing(model)) { stop("model is mandatory in clusterPredict.")}
  if(!is(model,"IClusterModel")) { stop("model must be an instance of IClusterModel or a derived class.")}
  
  # for data
  if(missing(data)) { stop("data is mandatory in clusterPredict.")}
  
  # cluster
  if (is(model, "ClusterMixedDataModel"))
  {
    nbComponent <- length(model@lcomponent)
    if(length(data) != nbComponent)
    { stop("data does not have the same number of component than the model")}
    # get nbsample of the first data set
    nbSample <- nrow(data[[1]])
    # check data
    for(i in 1:nbComponent)
    {
      if (nrow(data[[i]]) != nbSample)
      {stop("Error in clusterPredict validity. all data must have nbSample rows.")}
      if (.hasSlot(model@lcomponent[[i]],"plkj"))
      { 
        modelDim <- dim(model@lcomponent[[i]]@plkj)
        dim(model@lcomponent[[i]]@plkj) <- c(modelDim[1] * modelDim[2], modelDim[3])
      }
    }
    result = new("ClusterPredictMixedData", data, model, algo)
  }
  else
  {
    data = as.matrix(data)
    if (nrow(data)<1) { stop("data is empty in clusterPredict.")}
    # check dimension
    if (ncol(model@component@data) != ncol(data))
    { stop("data does not have the same number of variables than the model")}  
    result = new("ClusterPredict", data, model, algo)
    # transform plkj in matrix
    if ( is(model,"ClusterCategorical") )
    {
      modelDim <- dim(model@component@plkj)
      dim(model@component@plkj) <- c(modelDim[1] * modelDim[2], modelDim[3])
    } 
  }
  
  # start estimation of the models
  resFlag = .Call("clusterPredict", model, result, nbCore, PACKAGE="MixAll");
  if (resFlag != TRUE ) { cat("WARNING: An error occur during the predicting process");}
  
  # set plkj as array
  if ( is(model,"ClusterCategorical") )
  {
    nbVariable <- dim(model@component@plkj)[2]
    dim(model@component@plkj) <- c(model@component@nbModalities, model@nbCluster, nbVariable)
  }
  if (is(model, "ClusterMixedDataModel"))
  {
    nbComponent <- length(model@lcomponent)
    for(i in 1:nbComponent)
    {
      if (.hasSlot(model@lcomponent[[i]],"plkj"))
      {
        nbVariable <- dim(model@lcomponent[[i]]@plkj)[2]
        modelDim   <- c(model@lcomponent[[i]]@nbModalities, model@nbCluster, nbVariable)
        dim(model@lcomponent[[i]]@plkj) <- modelDim
      }
    }
  }
  #
  result
}


#' Class [\code{\linkS4class{ClusterPredict}}] for predicting 
#'
#' This class encapsulate the parameters for predicted data.
#'
#' @slot data  Matrix with the data set
#' @slot missing   Matrix with the indexes of the missing values
#'
#' @seealso [\code{\linkS4class{IClusterPredict}}] class
#' 
#' @examples
#'   getSlots("ClusterPredict")
#'
#' @author Serge Iovleff
#'
#' @name ClusterPredict
#' @rdname ClusterPredict-class
#' @aliases ClusterPredict-class
#' 
setClass(
  Class = "ClusterPredict",
  # members
  representation( data    = "matrix"
                , missing = "matrix"
                ),
  contains=c("IClusterPredict"),
  # validity function
  validity=function(object)
  {
    nbSample  <- object@nbSample
    nbCluster <- object@nbCluster
    
    # check data
    if (nrow(object@data) != nbSample)
    {stop("Error in ClusterPredict validity. data must have nbSample rows.")}
    return(TRUE)
  }
)

#' Initialize an instance of a MixAll S4 class.
#'
#' Initialization method of the [\code{\linkS4class{ClusterPredict}}] class.
#' Used internally in the 'MixAll' package.
#'
#' @rdname initialize-methods
#' @keywords internal
#'
setMethod(
  f="initialize",
  signature=c("ClusterPredict"),
  definition=function(.Object, data, model, algo = clusterAlgoPredict())
  {
    # for data
    if(missing(data)) { stop("data is mandatory in ClusterPredict.")}
    .Object@data     <- as.matrix(data)
    .Object@missing  <- which(is.na(.Object@data), arr.ind=TRUE);
    
    # for nbCluster
    if(missing(model)) { stop("model is mandatory in ClusterPredict.")}
    
    # create base class
    .Object <- callNextMethod(.Object, nrow(data), model, algo)
    
    # valid object
    validObject(.Object)
    
    # in the derived classes
    return(.Object)
  }
)

#' @rdname print-methods
#' @aliases print print,ClusterPredict-method
#'
setMethod(
  f="print",
  signature=c("ClusterPredict"),
  function(x,...)
  {
    cat("****************************************\n")
    if(length(x@data)!=0)
    {
      nrowShow <- min(10,nrow(x@data));
      ncolShow <- min(10,ncol(x@data));
      cat("* data (limited to 10 samples and 10 variables) =\n")
      print(format(x@data[1:nrowShow,1:ncolShow]),quote=FALSE)
    }
    callNextMethod()
    cat("****************************************\n")
  }
)

#' @rdname show-methods
#' @aliases show show,ClusterPredict-method
setMethod(
  f="show",
  signature=c("ClusterPredict"),
  function(object)
  {
    cat("****************************************\n")
    callNextMethod()
    cat("****************************************\n")
  }
)

#' @rdname summary-methods
#' @aliases summary summary,ClusterPredict-method
setMethod(
  f="summary",
  signature=c("ClusterPredict"),
  function(object,...)
  {
    cat("****************************************\n")
    callNextMethod()
    cat("****************************************\n")
  }
)

#' Class [\code{\linkS4class{ClusterPredictMixedData}}] for predicting 
#'
#' This class encapsulate the parameters for predicted data.
#'
#' @slot ldata    list of matrix with the data sets
#' @slot lmissing list of matrix with the indexes of the missing values
#'
#' @seealso [\code{\linkS4class{IClusterPredict}}] class
#' 
#' @examples
#'   getSlots("ClusterPredictMixedData")
#'
#' @author Serge Iovleff
#'
#' @name ClusterPredictMixedData
#' @rdname ClusterPredictMixedData-class
#' @aliases ClusterPredictMixedData-class
#' 
setClass(
  Class = "ClusterPredictMixedData",
  representation( ldata = "list", lmissing = "list"),
  contains=c("IClusterPredict"),
  # validity function
  validity=function(object)
  {
    nbSample   <- object@nbSample
    nbCluster  <- object@nbCluster
    nbComponent<- length(object@ldata)
    # check nbComponent
    if (nbComponent < 1)
    {stop("Error in ClusterPredictMixedData validity. There is no component")}
    # check data
    for(i in 1:nbComponent)
    {
      if (nrow(object@ldata[[i]]) != nbSample)
      {stop("Error in ClusterPredictMixedData validity. all data must have nbSample rows.")}
    }
    return(TRUE)
  }
)

#' Initialize an instance of a MixAll S4 class.
#'
#' Initialization method of the [\code{\linkS4class{ClusterPredictMixedData}}] class.
#' Used internally in the 'MixAll' package.
#'
#' @rdname initialize-methods
#' @keywords internal
#'
setMethod(
  f="initialize",
  signature=c("ClusterPredictMixedData"),
  definition=function(.Object, ldata, model, algo = clusterAlgoPredict())
  {
    # for ldata
    if(missing(ldata)) { stop("data is mandatory in ClusterPredictMixedData.")}
    .Object@ldata <- ldata

    # for model
    if(missing(model)) { stop("model is mandatory in ClusterPredictMixedData.")}

    # for nbComponent
    nbComponent<- length(ldata)
    if (nbComponent != length(model@lcomponent))
    { stop("Error in ClusterPredictMixedData. length(ldata) != length(model@lcomponent)")}
    
    # check nbSample using first data set
    nbSample <- nrow(ldata[[1]])
    .Object@lmissing <- vector("list", nbComponent)
    for(i in 1:nbComponent)
    {
      if (nrow(ldata[[i]]) != nbSample)
      { stop("Error in ClusterPredictMixedData initialize. all data must have the same number of rows.")}
      else
      { 
        if ( class(model@lcomponent[[i]]) == "ClusterCategoricalComponent")
        {
          data   <- as.data.frame(ldata[[i]])
          levels <- model@lcomponent[[i]]@levels
          for ( j in 1:length(data) )
          { 
            data[,j] <- as.integer(factor(data[,j]), levels = levels[[j]])}
            .Object@ldata[[i]] <- as.matrix(data)
        }
        else
        { 
          .Object@ldata[[i]] <- as.matrix(ldata[[i]])
        }
        .Object@lmissing[[i]] <- which(is.na(.Object@ldata[[i]]), arr.ind=TRUE)
      }
    }
    
    # create base class
    .Object <- callNextMethod(.Object, nbSample, model, algo)
    
    # valid object
    validObject(.Object)
    
    # in the derived classes
    return(.Object)
  }
)

#' @rdname print-methods
#' @aliases print print,ClusterPredictMixedData-method
#'
setMethod(
  f="print",
  signature=c("ClusterPredictMixedData"),
  function(x,...)
  {
    cat("****************************************\n")
    if(length(x@data)!=0)
    {
      nbComponent<- length(x@ldata)
      for(i in 1:nbComponent)
      {
        cat("* Component: ", i, "\n")
        nrowShow <- min(10,nrow(x@ldata[[i]]));
        ncolShow <- min(10,ncol(x@ldata[[i]]));
        cat("* data (limited to 10 samples and 10 variables) =\n")
        print(format(x@ldata[[i]][1:nrowShow,1:ncolShow]),quote=FALSE)
        cat("*\n")
      }
    }
    callNextMethod()
    cat("****************************************\n")
  }
)

#' @rdname show-methods
#' @aliases show show,ClusterPredictMixedData-method
setMethod(
  f="show",
  signature=c("ClusterPredictMixedData"),
  function(object)
  {
    cat("****************************************\n")
    callNextMethod()
    cat("****************************************\n")
  }
)

#' @rdname summary-methods
#' @aliases summary summary,ClusterPredictMixedData-method
setMethod(
  f="summary",
  signature=c("ClusterPredictMixedData"),
  function(object,...)
  {
    cat("****************************************\n")
    callNextMethod()
    cat("****************************************\n")
  }
)



