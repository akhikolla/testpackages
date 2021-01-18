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

#' Interface class [\code{\linkS4class{IClusterPredict}}] for predicting 
#'
#' Interface base class for predicting clusters
#'
#' @slot nbSample  Integer with the number of samples 
#' @slot nbCluster Integer with the number of cluster
#' @slot pk        Vector of size K with the proportions of each mixture.
#' @slot tik       Matrix of size \eqn{n \times K} with the posterior probability of
#' the ith individual to belong to kth cluster.
#' @slot lnFi      Vector of size n with the log-likelihood of the ith individuals.
#' @slot zi        Vector of integer of size n  with the attributed class label of the individuals
#' @slot algo      an instance of [\code{\linkS4class{ClusterAlgoPredict}}] 
#' @slot model     an instance of a (derived) [\code{\linkS4class{IClusterModel}}] 
#'
#'
#' @examples
#'   getSlots("IClusterPredict")
#'
#' @author Serge Iovleff
#'
#' @name IClusterPredict
#' @rdname IClusterPredict-class
#' @aliases IClusterPredict-class
#' 
setClass(
  Class = "IClusterPredict",
  # members
  representation( nbSample  = "numeric"
                , nbCluster = "numeric"
                , pk        = "numeric"
                , tik       = "matrix"
                , lnFi      = "numeric"
                , zi        = "integer"
                , algo      = "ClusterAlgoPredict"
                , model     = "IClusterModel"
                , "VIRTUAL"
                ),
  # validity function
  validity=function(object)
  {
    nbSample  <- object@nbSample
    nbCluster <- object@nbCluster
    
    # check nbSample
    if (round(nbSample) != nbSample)
    {stop("Error in IClusterPredict validity. nbSample must be an integer.")}
    
    # check nbCluster
    if (round(nbCluster)!= nbCluster)
    {stop("Error in IClusterPredict validity. nbCluster must be an integer.")}
    if( nbCluster < 1 )
    { stop("Error in IClusterPredict validity. nbCluster must be greater than 0.")}
    
    # check pk
    if (length(object@pk) != nbCluster)
    { stop("Error in IClusterPredict validity. pk must have length nbCluster.")}
    
    # check tik
    if (ncol(object@tik) != nbCluster)
    { stop("Error in IClusterPredict validity. tik must have nbCluster columns.")}
    if (nrow(object@tik) != nbSample)
    { stop("Error in IClusterPredict validity. tik must have nbSample rows.")}
    
    # check lnFi
    if (length(object@lnFi) != nbSample)
    { stop("Error in IClusterPredict validity. lnFi must have nbSample size.")}
    
    # check zi
    if (length(object@zi) != nbSample)
    { stop("Error in IClusterPredict validity. zi must have nbSample size.")}
    
    # check algo
    if (class(object@algo)[1] != "ClusterAlgoPredict")
    { stop("Error in IClusterPredict validity. algo must be an instance of ClusterPredictAlgo.")}
    return(TRUE)
  }
)

#' Initialize an instance of a MixAll S4 class.
#'
#' Initialization method of the [\code{\linkS4class{IClusterPredict}}] class.
#' Used internally in the 'MixAll' package.
#'
#' @rdname initialize-methods
#' @keywords internal
#'
setMethod(
  f="initialize",
  signature=c("IClusterPredict"),
  definition=function(.Object, nbSample, model, algo)
  {
    # for nbCluster
    if(missing(nbSample)) { stop("nbSample is mandatory in IClusterPredict.")}
    .Object@nbSample<-nbSample
    
    # for nbCluster
    if(missing(model)) { stop("model is mandatory in IClusterPredict.")}
    .Object@model<-model
    .Object@nbCluster<-model@nbCluster
    
    # for nbCluster
    if(missing(algo)) { stop("algo is mandatory in IClusterPredict.")}
    .Object@algo<-algo
    
    # create arrays
    .Object@pk   <- model@pk
    .Object@tik  <- matrix(1/.Object@nbCluster, .Object@nbSample, .Object@nbCluster)
    .Object@lnFi <- rep(0, .Object@nbSample)
    .Object@zi   <- as.integer(rep(1, .Object@nbSample))
    .Object@algo <- algo
    
    # valid object
    validObject(.Object)
    
    # in the derived classes
    return(.Object)
  }
)

#' @rdname print-methods
#' @aliases print print,IClusterPredict-method
#'
setMethod(
  f="print",
  signature=c("IClusterPredict"),
  function(x,...)
  {
    cat("* nbSample       = ", x@nbSample, "\n")
    cat("* nbCluster      = ", x@nbCluster, "\n")
    cat("* zi             =\n")
    print( format(x@zi), quote=FALSE)
  }
)

#' @rdname show-methods
#' @aliases show show,IClusterPredict-method
setMethod(
  f="show",
  signature=c("IClusterPredict"),
  function(object)
  {
    cat("* nbSample       = ", object@nbSample, "\n")
    cat("* nbCluster      = ", object@nbCluster, "\n")
    cat("* zi             =\n")
    print( format(object@zi), quote=FALSE)
  }
)

#' @rdname summary-methods
#' @aliases summary summary,IClusterPredict-method
setMethod(
  f="summary",
  signature=c("IClusterPredict"),
  function(object,...)
  {
    cat("* nbSample       = ", object@nbSample, "\n")
    cat("* nbCluster      = ", object@nbCluster, "\n")
  }
)


