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
NULL

#-----------------------------------------------------------------------
#' Create an instance of the [\code{\linkS4class{ClusterAlgoPredict}}] class
#'
#' A prediction algorithm is a two stage algorithm. In the first stage we perform
#' a Monte Carlo algorithm for simulating both missing values and latent class
#' variables. In the second stage, we simulate or impute missing values.
#' 
#' The epsilon value is not used when the algorithm is "SemiSEM".
#'
#' @param algo character string with the second stage estimation algorithm.
#' Possible values are "EM", "SemiSEM". Default value is "EM".
#' @param nbIterBurn Integer defining the maximal number of burning iterations.
#' Default value is 50.
#' @param nbIterLong Integer defining the maximal number of iterations.
#' Default value is 100.
#' @param epsilon Real defining the epsilon value for the algorithm. Not used
#' with "semiSEM" algorithms. Default value is 1.e-7.
#'
#' @examples
#' clusterAlgoPredict()
#' clusterAlgoPredict(algo="SemiSEM", nbIterBurn=0)
#' clusterAlgoPredict(algo="EM", epsilon = 1e-06)
#'
#' @return a [\code{\linkS4class{ClusterAlgoPredict}}] object
#' @author Serge Iovleff
#'
#'
clusterAlgoPredict <- function( algo="EM", nbIterBurn = 50, nbIterLong = 100, epsilon=1e-07)
{
  # for algo
  if (!is.character(algo) )
  { stop("algo must be a string.")}
  algo = toupper(algo)
  
  # check algo
  if ( sum(algo %in% c("EM","SEMISEM")) != 1 )
  {  stop("algo is not valid. See ?clusterAlgo for the list of available algorithms.")}

  # check nbIterBurn
  if (!is.numeric(nbIterBurn))
  {stop("nbIterBurn must be an integer.")}
  if (round(nbIterBurn)!= nbIterBurn)
  {stop("nbIterBurn must be an integer.")}
  if( nbIterBurn < 0 ) # can be zero (no iterations)
  { stop("nbIterBurn must be positive or zero.")}
  
  # check nbIterLong
  if (!is.numeric(nbIterLong))
  {stop("nbIterLong must be an integer.")}
  if (round(nbIterLong)!= nbIterLong)
  {stop("nbIterLong must be an integer.")}
  if( nbIterLong < 0 ) # can be zero (no iterations)
  { stop("nbIterLong must be positive or zero.")}
  
  # check epsilon
  if (algo == "EM")
  {
    if (!is.double(epsilon) )
    {  stop("epsilon must be a scalar.")}
    if( epsilon < 0.)
    {  stop("epsilon must be positive.")}
  }
  new("ClusterAlgoPredict", algo=algo, nbIterBurn=nbIterBurn, nbIterLong=nbIterLong, epsilon=epsilon)
}

#-----------------------------------------------------------------------
#' @title [\code{\linkS4class{ClusterAlgoPredict}}] class for predict algorithm.
#'
#' @description This class encapsulates the parameters of prediction methods.
#'
#' @slot algo A character string with the algorithm.
#' Possible values: "EM", "SemiSEM". Default value: "SemiSEM".
#' @slot nbIterBurn Integer defining the number of burning iterations. Default value is 50.
#' @slot nbIterLong Integer defining the number of iterations. Default value is 100.
#' @slot epsilon real defining the epsilon value for the long algorithm. epsilon is
#' note used if \code{algo} is "SemiSEM". Default value: 1e-07.
#'
#' @examples
#' getSlots("ClusterAlgoPredict")
#' new("ClusterAlgoPredict")
#' new("ClusterAlgoPredict", algo="SemiSEM", nbIterBurn=10)
#'
#' @name ClusterAlgoPredict
#' @rdname ClusterAlgoPredict-class
#' @aliases ClusterAlgoPredict-class
#' 
setClass (
  Class= "ClusterAlgoPredict",
  representation(algo="character"
                ,nbIterBurn = "numeric"
                ,nbIterLong = "numeric"
                ,epsilon = "numeric"
                ),
  prototype=list(algo="SemiSEM", nbIterBurn = 50, nbIterLong = 100, epsilon = 1e-07),
  # validity function
  validity = function(object)
  {
    # for algo
    if (!is.character(object@algo) )
    { stop("algo must be a string.")}
    object@algo = toupper(object@algo)
    if ( sum(object@algo %in% c("EM","SEMISEM")) != 1 )
    { stop("Algorithm is not valid. See ?ClusterAlgoPredict for the list of available algorithms.")}

    # for nbIterBurn
    if (!is.numeric(object@nbIterBurn))
    { stop("nbIterBurn must be an integer.")}
    if (round(object@nbIterBurn)!=object@nbIterBurn)
    { stop("nbIterBurn must be an integer.")}
    if( object@nbIterBurn < 0 ) # can be zero (no iterations)
    { stop("nbIterBurn must be positive or zero.")}

    # for nbIterLong
    if (!is.numeric(object@nbIterLong))
    { stop("nbIterLong must be an integer.")}
    if (round(object@nbIterLong)!=object@nbIterLong)
    { stop("nbIterLong must be an integer.")}
    if( object@nbIterLong <0 ) # can be zero (no iterations)
    { stop("nbIterLong must be positive or zero.")}
    
    # for epsilon
    if (object@algo != "SemiSEM")
    {
      if (!is.double(object@epsilon) )
      { stop("epsilon must be a scalar.")}
      if( object@epsilon < 0.)
      { stop("epsilon must be positive.")}
    }
    return(TRUE)
  }
)

#' Initialize an instance of a MixAll S4 class.
#'
#' Initialization method of the [\code{\linkS4class{ClusterAlgoPredict}}] class.
#' Used internally in the `MixAll' package.
#'
#' @keywords internal
#' @rdname initialize-methods
setMethod(
    f="initialize",
    signature=c("ClusterAlgoPredict"),
    definition=function(.Object, algo, nbIterBurn, nbIterLong, epsilon)
    {
      # for algo
      if(missing(algo)) {.Object@algo<-"SemiSEM"}
      else  {.Object@algo<-algo}
      
      # for nbIterBurn
      if(missing(nbIterBurn)){ .Object@nbIterBurn<-50 }
      else{.Object@nbIterBurn<-nbIterBurn}
      
      # for nbIterLong
      if(missing(nbIterLong)){ .Object@nbIterLong<-50 }
      else{.Object@nbIterLong<-nbIterLong}
      
      # for epsilon
      if( missing(epsilon) ){ .Object@epsilon<-1e-07 }
      else{.Object@epsilon<-epsilon}
      # validate
      validObject(.Object)
      return(.Object)
    }
)

#' @rdname print-methods
#' @aliases print print-algo,ClusterAlgoPredict-method
setMethod(
    f="print",
    signature=c("ClusterAlgoPredict"),
    function(x,...){
      cat("****************************************\n")
      cat("*** MixAll ClusterAlgoPredict:\n")
      cat("* algorithm       = ", x@algo, "\n")
      cat("* burn iterations = ", x@nbIterBurn, "\n")
      cat("* long iterations = ", x@nbIterLong, "\n")
      cat("* epsilon         = ", x@epsilon, "\n")
      cat("****************************************\n")
    }
)

#' @rdname show-methods
#' @aliases show show-algo,ClusterAlgoPredict-method
setMethod(
    f="show",
    signature=c("ClusterAlgoPredict"),
    function(object){
      cat("****************************************\n")
      cat("*** MixAll ClusterAlgoPredict:\n")
      cat("* algorithm       = ", object@algo, "\n")
      cat("* burn iterations = ", object@nbIterBurn, "\n")
      cat("* long iterations = ", object@nbIterLong, "\n")
      cat("* epsilon         = ", object@epsilon, "\n")
      cat("****************************************\n")
    }
)

#' @rdname extract-methods
#' @aliases [,ClusterAlgoPredict-method
setMethod(
  f="[",
  signature(x = "ClusterAlgoPredict"),
  definition=function(x, i, j, drop){
    if ( missing(j) ){
      switch(EXPR=i,
          "algo"={return(x@algo)},
          "nbIterBurn"={return(x@nbIterBurn)},
          "nbIterLong"={return(x@nbIterLong)},
          "epsilon"={return(x@epsilon)},
          stop("This attribute doesn't exist !")
      )
    }else{
      stop("This attribute is not a list !")
    }
  }
)

#' @name [
#' @rdname extract-methods
#' @aliases [<-,ClusterAlgoPredict-method
setReplaceMethod(
    f="[",
    signature(x = "ClusterAlgoPredict"),
    definition=function(x,i,j,value){
      if ( missing(j) )
      {
        switch(EXPR=i,
            "algo"={x@algo<-value},
            "nbIterBurn"={x@nbIterBurn<-value},
            "nbIterLong"={x@nbIterLong<-value},
            "epsilon"={x@epsilon<-value},
            stop("This attribute doesn't exist !")
        )
      }else{
        stop("This attribute is not a list !")
      }
      validObject(x)
      return(x)
    }
)
