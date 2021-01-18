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
#' @include kmmNames.R IClusterModel.R
NULL

#-----------------------------------------------------------------------
#' Create an instance of [\code{\linkS4class{ClusterStrategy}}] class
#'
#' A strategy is a multistage empirical process for finding a
#' good estimate in the clustering estimation process.
#'
#' A strategy is a way to find a good estimate of the parameters of a kernel
#' mixture model when using an EM algorithm or its variants. A ``try'' of
#' kmmStrategy is composed of three stages
#' \itemize{
#'   \item \code{nbShortRun} short iterations of the initialization step and
#'    of the \code{EM}, \code{CEM} or \code{SEM} algorithm.
#'   \item \code{nbInit} initializations using the [\code{\link{clusterInit}}]
#'   method.
#'   \item A long run of the \code{EM}, \code{CEM} or \code{SEM} algorithm.
#' }
#' For example if \code{nbInit} is 5 and \code{nbShortRun} is also 5, there will
#' be 5 times 5 models initialized. Five time, the best model (in the likelihood sense)
#' will be ameliorated using a short run. Among the 5 models ameliorated one will be
#' estimated until convergence using a long run. In total there is 25 initializations.
#'
#' The whole process can be repeated at least \code{nbTry} times. If a try
#' success, the estimated model is returned, otherwise an empty model is returned.
#'
#' @param nbTry Integer defining the number of estimation to attempt.
#'
#' @param nbInit Integer defining the number of initialization to try. Default value: 3.
#' @param initMethod Character string with the initialization method, see [\code{\link{clusterInit}}]$
#' for possible values. Default value: "class".
#' @param initAlgo Character string with the algorithm to use in the initialization stage,
#' [\code{\link{clusterAlgo}}] for possible values. Default value: "EM".
#' @param nbInitIteration Integer defining the maximal number of iterations in
#' initialization algorithm if \code{initAlgo} = "EM" or "CEM", the number of iterations
#' if \code{initAlgo} = "SEM". Default value: 20.
#' @param initEpsilon Real defining the epsilon value for the initialization algorithm.
#' Not used if  \code{initAlgo} = "SEM". Default value: 0.01.
#'
#' @param nbShortRun Integer defining the number of short run to try
#' (the strategy launch an initialization before each short run). Default value: 5.
#' @param shortRunAlgo A character string with the algorithm to use in the short run stage.
#' Default value: "EM".
#' @param nbShortIteration Integer defining the maximal number of iterations during
#' sa hort run if \code{shortRunAlgo} = "EM" or "CEM", the number of iterations
#' if \code{shortRunAlgo} = "SEM". Default value: 100.
#' @param shortEpsilon Real defining the epsilon value for the algorithm. Not used
#' if \code{shortRunAlgo} = "SEM". Default value: 1e-04.
#'
#' @param longRunAlgo A character string with the algorithm to use in the long run stage.
#' Default value: "EM".
#' @param nbLongIteration Integer defining the maximal number of iterations during
#' a long run algorithm if \code{longRunAlgo} = "EM" or "CEM", the number of iterations
#' if \code{longRunAlgo} = "SEM". Default value: 1000.
#' @param longEpsilon Real defining the epsilon value for the algorithm.
#' Nor used if \code{longRunAlgo} = "SEM". Default value: 1e-07.
#'
#' @examples
#'    kmmStrategy()
#'    kmmStrategy(longRunAlgo= "CEM", nbLongIteration=100)
#'    kmmStrategy(nbTry = 1, nbInit= 1, shortRunAlgo= "EM", nbShortIteration=100)
#'
#' @return a [\code{\linkS4class{ClusterStrategy}}] object
#' @author Serge Iovleff
#'
kmmStrategy <- function( nbTry =1
                       , nbInit= 5, initMethod="class", initAlgo= "EM", nbInitIteration=20, initEpsilon=0.01
                       , nbShortRun= 5, shortRunAlgo= "EM", nbShortIteration=100, shortEpsilon=1e-04
                       , longRunAlgo= "EM", nbLongIteration=1000, longEpsilon=1e-07
                       )
{
  # check initAlgo
  initAlgo = toupper(initAlgo)
  if ( sum(initAlgo %in% c("EM","SEM","CEM")) != 1 )
  { stop("initAlgo must be EM, CEM or SEM")}
  # check shortRunAlgo
  shortRunAlgo = toupper(shortRunAlgo)
  if ( sum(initAlgo %in% c("EM","SEM","CEM")) != 1 )
  { stop("shortRunAlgo must be EM, CEM or SEM")}
  # check longRunAlgo
  longRunAlgo = toupper(longRunAlgo)
  if ( sum(longRunAlgo %in% c("EM","SEM","CEM")) != 1 )
  { stop("longRunAlgo must be EM, CEM or SEM")}
  # create init
  init = clusterInit(initMethod, nbInit, initAlgo, nbInitIteration, initEpsilon);
  # create shortAlgo
  shortAlgo = clusterAlgo(shortRunAlgo, nbShortIteration, shortEpsilon);
  # create longAlgo
  longAlgo = clusterAlgo(longRunAlgo, nbLongIteration, longEpsilon);
  # create strategy
  new("ClusterStrategy", nbTry =nbTry, nbShortRun =nbShortRun, initMethod =init, shortAlgo =shortAlgo, longAlgo =longAlgo);
}


#-----------------------------------------------------------------------
#' Create an instance of the [\code{\linkS4class{KmmModel}}] class
#'
#' This function computes the optimal kernel mixture model (KMM) according
#' to the [\code{criterion}] among the number of clusters given in
#' [\code{nbCluster}], using the strategy specified in [\code{strategy}].
#'
#' @param data frame or matrix containing the data. Rows correspond to observations
#' and columns correspond to variables.
#' @param dim integer giving the dimension of the Gaussian density. Default is 10.
#' @param nbCluster  [\code{\link{vector}}] listing the number of clusters to test.
#' @param models [\code{\link{vector}}] of model names to run. By default only
#' "kmm_pk_s" is estimated. All the model names are given by the method
#' [\code{\link{kmmNames}}].
#' @param kernelName string with a kernel name. Possible values:
#' "Gaussian", "polynomial", "Laplace", "linear", "rationalQuadratic_", "Hamming".
#' Default is "Gaussian".
#' @param kernelParameters [\code{\link{vector}}] with the parameters of
#' the chosen kernel. Default is c(1).
#' @param kernelComputation [\code{\link{logical}}] parameter. Should be \code{TRUE}
#' if the Gram matrix is to be computed (faster but can be memory consuming), \code{FALSE}
#' otherwise (times consuming). Default is \code{TRUE}. Recall that Gram matrix
#' is a square matrix of size nbSample.
#' @param strategy a [\code{\linkS4class{ClusterStrategy}}] object containing
#' the strategy to run. [\code{\link{kmmStrategy}}]() method by default.
#' @param criterion character defining the criterion to select the best model.
#' The best model is the one with the lowest criterion value.
#' Possible values: "BIC", "AIC", "ICL", "ML". Default is "ICL".
#' @param nbCore integer defining the number of processor to use (default is 1, 0 for all).
#'
#' @note in KmmModel instance returned, the gram matrix is computed if and only
#' if kernelComputation is \code{TRUE}.
#' 
#' @examples
#' ## A quantitative example with the famous bulls eye model
#' data(bullsEye)
#' ## estimate model
#' model <- kmm( data=bullsEye, nbCluster=2:3, models= "kmm_pk_s")
#'
#'
#' ## get summary
#' summary(model)
#' ## use graphics functions
#' \dontrun{
#' plot(model)
#' }
#'
#' @return An instance of the [\code{\linkS4class{KmmModel}}] class.
#' @author Serge Iovleff
#'
#'
kmm <- function( data, nbCluster=2
               , dim = 10, models = "kmm_pk_s"
               , kernelName = "Gaussian", kernelParameters = c(1), kernelComputation = TRUE
               , strategy=kmmStrategy()
               , criterion="ICL"
               , nbCore = 1)
{
  # check data
  if (missing(data)) { stop("data is mandatory in kmm")}
  data <- as.matrix(data)
  if (ncol(data) < 1) { stop("Error: empty data set")}
  
  # check nbCluster
  if (length(nbCluster) == 0) { stop("You must give at least one number of clusters")}
  nbClusterMin = min(nbCluster);
  nbClusterMax = max(nbCluster);
  if (nbClusterMin < 1) { stop("The number of clusters must be greater or equal to 1")}  
  if (nrow(data) <= 3*nbClusterMax)
  { stop("There is too much clusters or not enough individuals in the data set")}
  
  # check dim
  if (dim < 1) { stop("dim must be greater or equal to 1")}
  
  # check criterion
  criterion = toupper(criterion)
  if(sum(criterion %in% c("BIC","AIC", "ICL", "ML")) != 1)
  { stop("criterion is not valid. See ?kmm for the list of valid criterion")}
  
  # check model names
  if (!kmmValidModelNames(models))
  { stop("models is not valid. See ?kmmNames for the list of valid model names")}
  
  # check kernelName
  if (!kmmValidKernelNames(kernelName))
  { stop("kernelName is not valid. See ?kmm for the list of valid kernel name")}
  if (is.null(kernelParameters)) { kernelParameters = c(1)}
 
  # check kernelComputation
  if (!is.logical(kernelComputation))
  { stop("kernelComputation is not boolean")}
  
  # check strategy
  if(class(strategy)[1] != "ClusterStrategy")
  {stop("strategy is not a ClusterStrategy class (must be an instance of the class ClusterStrategy)")}
  validObject(strategy);

  # Create model
  model = new( "KmmModel", data
             , dim=dim
             , kernelName = kernelName, kernelParameters = kernelParameters, kernelComputation=kernelComputation)
  model@strategy = strategy;
  model@criterionName = criterion;
  
  # start estimation of the models
  resFlag <- FALSE
  resFlag = .Call("kmm", model, nbCluster, models, nbCore, PACKAGE="MixAll")
  if (resFlag != TRUE )
  { 
    msg_error <- model@msg_error
    cat("WARNING: An error occur during the clustering process.\n")
    cat("Collected error:",msg_error)
  }
  model
}

#-----------------------------------------------------------------------
#' Definition of the [\code{\linkS4class{KmmComponent}}] class
#'
#' This class defines a kernel component of a mixture Model. It inherits
#' from [\code{\linkS4class{IClusterComponent}}].
#'
#' @slot dim    Vector with the dimension of the kth cluster
#' @slot sigma2 Vector with the standard deviation in the kth cluster.
#' @slot gram   Matrix storing the gram matrix if its computation is needed
#' @slot kernelName string with the name of the kernel to use. Possible values:
#' "Gaussian", "polynomial", "Laplace", "linear","rationalQuadratic", "Hamming".
#' Default is "Gaussian".
#' @slot kernelParameters vector with the parameters of the kernel.
#' @slot kernelComputation boolean value set as \code{TRUE} if Gram matrix is to be computed
#' \code{FALSE} othewise. Default is \code{TRUE}.
#'
#' @seealso [\code{\linkS4class{IClusterComponent}}] class
#'
#' @examples
#' getSlots("KmmComponent")
#'
#' @author Serge Iovleff
#'
#' @name KmmComponent
#' @rdname KmmComponent-class
#' @aliases KmmComponent-class
#'
setClass(
  Class = "KmmComponent",
  representation( sigma2 = "vector"
                , dim = "vector"
                , kernelName = "character"
                , kernelParameters = "vector"
                , kernelComputation = "logical"
                , gram = "matrix"),
  contains=c("IClusterComponent"),
  validity=function(object)
  {
    if ( length(object@sigma2) == 0) { stop("sigma2 must be a vector of length > 0.")}
    if ( length(object@dim)    == 0) { stop("dim must be a vector of length > 0")}
    
    if (!kmmValidModelNames(object@modelName))
    {stop(paste("Invalid model name: ", object@modelName,". See ?KmmModel for the list of valid model name",sep=""))}
    if (!kmmValidKernelNames(object@kernelName))
    { stop(paste("kernelName is not valid: ",object@kernelName,". See ?KmmModel for the list of valid kernel name",sep=""))}
    
    return(TRUE)
  }
)
#' Initialize an instance of a KmmComponent S4 class.
#'
#' Initialization method of the [\code{\linkS4class{KmmComponent}}] class.
#' Used internally in the 'MixAll' package.
#'
#' @rdname initialize-methods
#' @keywords internal
setMethod(
    f="initialize",
    signature=c("KmmComponent"),
    definition=function( .Object, data, dim =10
                       , nbCluster=2
                       , modelName="kmm_pk_s"
                       , kernelName = "Gaussian"
                       , kernelParameters = c(1)
                       , kernelComputation=TRUE
                       )
    {
      # check data
      if(missing(data)) { stop("data is mandatory in KmmComponent"); }
      
      # check dim
      if(is.null(dim)) { dim = 10}
      else
      {
        if (length(dim) != 1)
        { stop("dim must be of length 1")}
      }
      .Object@dim = rep(dim, nbCluster)
      
      # check model name
      if(is.null(modelName)) { modelName="kmm_pk_s"}
      else
      {
        if (length(modelName) != 1)
        { stop("modelName must be of length 1.")}
        if(!kmmValidModelNames(modelName))
        { stop("modelName is invalid. See ?kmmNames for the list of valid model name");}
      }
      
      # check kernelName
      if(is.null(kernelName)) { kernelName="Gaussian"}
      else
      {
        if (length(kernelName) != 1)
        { stop("kernelName must be of length 1.")}
        if (!kmmValidKernelNames(kernelName))
        { stop(paste("kernelName is not valid: ",kernelName,". See ?KmmNames for the list of valid kernel name",sep=""))}
      }
      .Object@kernelName <- kernelName
      
      # check kernelParameters
      if (is.null(kernelParameters)) { kernelParameters = c(1);}
      else
      {
        if (length(kernelParameters) < 1)
        { stop("kernelParameters must be at least of length 1.")}
      }
      .Object@kernelParameters <- kernelParameters;
      
      # check kernelComputation
      if(is.null(kernelComputation)) { kernelComputation=TRUE}
      else
      {
        if (!is.logical(kernelComputation))
        { stop("kernelComputation must be boolean");}
      }
      .Object@kernelComputation <- kernelComputation;
      
      # create slots
      .Object@sigma2   = rep(1., nbCluster)
      .Object@gram     = matrix(nrow=0, ncol=0)
      
      # call base class initialize (use data as data)
      .Object <- callNextMethod(.Object, data, modelName)
      # validate
      validObject(.Object)
      return(.Object)
    }
)

#' @rdname extract-methods
#' @aliases [,KmmComponent-method
setMethod(
    f="[",
    signature(x = "KmmComponent"),
    definition=function(x, i, j, drop)
    {
      if ( missing(j) )
      {
        switch(EXPR=i,
            "sigma2" = {return(x@sigma2)},
            "dim"    = {return(x@dim)},
            stop("This attribute doesn't exist !")
        )
      }
      else
      {
        if (!is.numeric(j)) {stop("j must be an integer.")}
        if (round(j)!=j)    {stop("j must be an integer.")}
        switch(EXPR=i,
            "sigma2"  = {return(x@sigma2[j,])},
            "dim"    = {return(x@dim[j,])},
            stop("This attribute doesn't exist !")
        )
      }
    }
)

#' @rdname print-methods
#' @aliases print print,KmmComponent-method
#'
setMethod(
  signature=c("KmmComponent"),
  f="print",
  function(x,k,...)
  {
    cat("* sigma2  = ", format(x@sigma2[k]), "\n")
    cat("* dim    = ", format(x@dim[k]), "\n")
  }
)

#' @rdname show-methods
#' @aliases show-KmmComponent,KmmComponent,KmmComponent-method
setMethod(
  f="show",
  signature=c("KmmComponent"),
  function(object)
  {
    cat("* sigma2  = ", format(object@sigma2), "\n")
    cat("* dim    = " , format(object@dim), "\n")
  }
)

#-----------------------------------------------------------------------
#' Definition of the [\code{\linkS4class{KmmModel}}] class
#'
#' This class defines a Kernel mixture Model (KMM).
#'
#' This class inherits from the [\code{\linkS4class{IClusterModel}}] virtual class.
#' A KMM is a mixture model of the form:
#' \deqn{
#'   f({x}|\boldsymbol{\theta})
#'   =\sum_{k=1}^K p_k \prod_{j=1}^d \phi(x_j;\sigma^2_{k})
#'    \quad x \in {R}^d.
#' }
#' Some constraints can be added to the variances in order to reduce the number
#' of parameters.
#'
#' @slot component  A [\code{\linkS4class{KmmComponent}}] with the
#' dimension and standard deviation of the kernel mixture model.
#' @seealso [\code{\linkS4class{IClusterModel}}] class
#'
#' @examples
#' getSlots("KmmModel")
#' data(bullsEye)
#' new("KmmModel", data=bullsEye)
#'
#' @author Serge Iovleff
#'
#' @name KmmModel
#' @rdname KmmModel-class
#' @aliases KmmModel-class
#'
setClass(
  Class = "KmmModel",
  representation( component = "KmmComponent"),
  contains=c("IClusterModel"),
  validity=function(object)
  {
    if (length(object@component@dim)!=object@nbCluster)
    {stop("dim must have nbCluster length.")}
    if (length(object@component@sigma2)!=object@nbCluster)
    {stop("sigma2 must have nbCluster length.")}
    
    if (!kmmValidModelNames(object@component@modelName))
    {stop(paste("Invalid kernel mixture model name:", object@component@modelName,". See ?KmmModel for the list of valid kernel name",sep=""))}
    
    # check kernelName
    if (length(object@component@kernelName) != 1)
    { stop("kernelName must be of length 1. See ?KmmModel for the list of valid kernel name")}
    if (!kmmValidKernelNames(object@component@kernelName))
    { stop(paste("kernelName is not valid",object@component@kernelName,". See ?KmmModel for the list of valid kernel name",sep=""))}
    return(TRUE)
  }
)

#' Initialize an instance of a MixAll S4 class.
#'
#' Initialization method of the [\code{\linkS4class{KmmModel}}] class.
#' Used internally in the 'MixAll' package.
#'
#' @rdname initialize-methods
#' @keywords internal
setMethod(
    f="initialize",
    signature=c("KmmModel"),
    definition=function(.Object, data, nbCluster=2
                       , modelName = "kmm_pk_s", dim= 10
                       , kernelName = "Gaussian"
                       , kernelParameters = c(1)
                       , kernelComputation=TRUE
                       )
    {
      # check data
      if(missing(data)) {stop("data is mandatory in KmmModel.")}
      
      # check model name
      if(is.null(modelName)) { modelName="kmm_pk_s";}
      else
      { if(!kmmValidModelNames(modelName)) { stop("modelName is invalid");} }
            
      # initialize component
      .Object@component = new( "KmmComponent", data, dim, nbCluster
                             , modelName
                             , kernelName, kernelParameters, kernelComputation);
      .Object <- callNextMethod(.Object, nrow(data), nbCluster);
      # validate
      validObject(.Object);
      return(.Object);
    }
)

#' @rdname print-methods
#' @aliases print print,KmmModel-method
#'
setMethod(
  f="print",
  signature=c("KmmModel"),
  function(x,...){
    cat("****************************************\n")
    callNextMethod();
    cat("* KernelName = ", x@component@kernelName, "\n")
    cat("* KernelParameters = ", format(x@component@kernelParameters), "\n")
    cat("****************************************\n")
    for(k in 1:length(x@pk))
    {
      cat("*** Cluster: ",k,"\n")
      cat("* Proportion = ", format(x@pk[k]), "\n")
      print(x@component, k);
      cat("****************************************\n")
    }
  }
)

#' @rdname show-methods
#' @aliases show-KmmModel,KmmModel,KmmModel-method
setMethod(
  f="show",
  signature=c("KmmModel"),
  function(object)
  {
    cat("****************************************\n")
    callNextMethod();
    show(object@component);
    cat("* KernelName = ", object@component@kernelName, "\n")
    cat("* KernelParameters = ", format(object@component@kernelParameters), "\n")
    cat("****************************************\n")
    for(k in 1:length(object@pk))
    {
      cat("*** Cluster: ",k,"\n")
      cat("* Proportion = ", format(object@pk[k]), "\n")
      print(object@component, k);
      cat("****************************************\n")
    }
  }
)

#' @rdname summary-methods
#' @aliases summary summary,KmmModel-method
#'
setMethod(
  f="summary",
  signature=c("KmmModel"),
  function(object, ...)
  {
    cat("****************************************\n")
    callNextMethod()
    cat("* KernelName     = ", object@component@kernelName, "\n")
    summary(object@component);
    cat("****************************************\n")
  }
)

#' Plotting of a class [\code{\linkS4class{KmmModel}}]
#'
#' Plotting data from a [\code{\linkS4class{KmmModel}}] object
#' using the estimated parameters and partition.
#'
#' @param x an object of class [\code{\linkS4class{KmmModel}}]
#' @param y a list of variables to plot (subset). Variables names or indices.
#' If missing all the variables are represented.
#' @param ... further arguments passed to or from other methods
#'
#' @aliases plot-KmmModel
#' @docType methods
#' @rdname plot-KmmModel-method
#'
#'
#' @seealso \code{\link{plot}}
#' @examples
#' \dontrun{
#'  ## the bull eyes data set
#'   data(bullsEye)
#'   model <- kmm( bullsEye, 2, models= "kmm_pk_s")
#'   plot(model)
#'   }
#'
setMethod(
    f="plot",
    signature=c("KmmModel"),
    function(x, y, ...)
    {
      # total number of variable in the data set
      data = x@component@data
      nbVariable = ncol(data);
      # no y => display all variables
      if (missing(y)) { y=1:nbVariable; }
      else # perform some check
      {
        if (is.numeric(y)) # numbers of the columns to plot are given
        {
          if (max(y)>nbVariable)
            stop("In plot, y indices mismatch the data dimension")
        }
        else # names of the variables to plot are given
        {
          if ( sum(y %in% colnames(data)) != length(y) )
          { stop(cat("In plot, unknown variables: ", paste(y[which(!(y %in% colnames(data)))]),"\n"))}
        }
      }
      # scatter plot
      plot(as.data.frame(data[,y]), col = x@zi+2)
    }
)

#' Plotting of a class [\code{\linkS4class{KmmComponent}}]
#'
#' Plotting data from a [\code{\linkS4class{KmmComponent}}] object
#' using the estimated partition.
#'
#' @param x an object of class [\code{\linkS4class{KmmComponent}}]
#' @param y a vector with partitions
#' @param ... further arguments passed to or from other methods
#'
#' @aliases plot-KmmComponent
#' @docType methods
#' @rdname plot-KmmComponent-method
#'
#'
#' @seealso \code{\link{plot}}
#' @examples
#' \dontrun{
#'  ## the bull eyes data set
#'   data(bullsEye)
#'   model <- kmm( bullsEye, 2, models= "kmm_pk_s")
#'   plot(model)
#'   }
#'
setMethod(
  f="plot",
  signature=c("KmmComponent"),
  function(x, y, ...)
  {
    # total number of variable in the data set
    data = x@data
    nbVariable = ncol(data);
    # no y => no class
    if (missing(y)) { y=rep(1,nrow(data)) }
    # scatter plot using data frame plot
    plot(as.data.frame(data), col = y+2)
  }
)

