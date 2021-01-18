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
#' Create an instance of the [\code{\linkS4class{LearnAlgo}}] class
#'
#' There is two algorithms and two stopping rules possibles for a learning
#' algorithm.
#' \itemize{
#'        \item Algorithms:
#'           \itemize{
#'               \item \code{Impute} {Impute the missing values during the iterations}
#'               \item \code{Simul} {Simulate the missing values during the iterations}
#'           }
#'        \item Stopping rules:
#'           \itemize{
#'               \item \code{nbIteration} {Set the maximum number of iterations.}
#'               \item \code{epsilon} {Set relative increase of the log-likelihood criterion.}
#'           }
#'        \item Default values are \eqn{200} \code{nbIteration} of \code{Simul}.
#' }
#' The \code{epsilon} value is not used when the algorithm is "Simul". It is worth noting
#' that if there is no missing values, the method should be "Impute" and nbIteration
#' should be set to 1!
#'
#' @param algo character string with the estimation algorithm.
#' Possible values are "Simul", "Impute". Default value is "Simul".
#' @param nbIteration Integer defining the maximal number of iterations. Default value is 200.
#' @param epsilon Real defining the epsilon value for the algorithm. Not used
#'  by the "Simul" algorithm. Default value is 1.e-7.
#'
#' @examples
#' learnAlgo()
#' learnAlgo(algo="simul", nbIteration=50)
#' learnAlgo(algo="impute", epsilon = 1e-06)
#'
#' @return a [\code{\linkS4class{LearnAlgo}}] object
#' @author Serge Iovleff
#'
#'
learnAlgo <- function( algo="Simul", nbIteration=200, epsilon=1e-07)
{
  # check algo
  if (!is.character(algo) )
  { stop("algo must be a string.")}
  algo = toupper(algo)
  if ( sum(algo %in% c("SIMUL","IMPUTE")) != 1 )
  { stop("algo is not valid. See ?learnAlgo for the list of available algorithms.")}
  # check nbIteration
  if (!is.numeric(nbIteration))
  {stop("nbIteration must be an integer.")}
  if (round(nbIteration)!= nbIteration)
  {stop("nbIteration must be an integer.")}
  if( nbIteration < 0 ) # can be zero (no iterations)
  { stop("nbIteration must be positive or zero.")}
  # check epsilon
  if (algo != "SIMUL")
  {
    if (!is.double(epsilon) )
    {  stop("epsilon must be a scalar.")}
    if( epsilon < 0.)
    {  stop("epsilon must be positive.")}
  }
  new("LearnAlgo", algo=algo, nbIteration=nbIteration, epsilon=epsilon)
}

#-----------------------------------------------------------------------
#' @title [\code{\linkS4class{LearnAlgo}}] class for Cluster algorithms.
#'
#' @description
#' This class encapsulates the parameters of clustering estimation algorithms
#' methods.
#'
#' @slot algo A character string with the algorithm.
#' Possible values: "Simul", "Impute. Default value: "Simul".
#' @slot nbIteration Integer defining the maximal number of iterations. Default value: 200.
#' @slot epsilon real defining the epsilon value for the algorithm. epsilon is
#' note used if \code{algo} is "Simul". Default value: 1e-07.
#'
#' @examples
#' getSlots("LearnAlgo")
#' new("LearnAlgo")
#' new("LearnAlgo", algo="Impute", nbIteration=100)
#'
#' @name LearnAlgo
#' @rdname LearnAlgo-class
#' @aliases LearnAlgo-class
#'
setClass (
  Class= "LearnAlgo",
  representation(algo="character", nbIteration = "numeric", epsilon = "numeric"),
  prototype=list(algo="Simul", nbIteration = 200, epsilon = 1e-07),
  # validity function
  validity = function(object)
  {
    # for algo
    if (!is.character(object@algo) )
    { stop("algo must be a string.")}
    object@algo = toupper(object@algo)
    # for algo
    if ( sum(object@algo %in% c("SIMUL","IMPUTE")) != 1 )
    {  stop("Algorithm is not valid. See ?LearnAlgo for the list of available algorithms.")}
    # for nbIteration
    if (!is.numeric(object@nbIteration))
    {stop("nbIteration must be an integer.")}
    if (round(object@nbIteration)!= object@nbIteration)
    {stop("nbIteration must be an integer.")}
    if( object@nbIteration < 0 ) # can be zero (no iterations)
    {  stop("nbIteration must be positive or zero.")}
    # for epsilon
    if ( object@algo != "Simul")
    {
      if (!is.double(object@epsilon) )
      {  stop("epsilon must be a scalar.")}
      if( object@epsilon < 0.)
      {  stop("epsilon must be positive.")}
    }
    return(TRUE)
  }
)

#' Initialize an instance of a MixAll S4 class.
#'
#' Initialization method of the [\code{\linkS4class{LearnAlgo}}] class.
#' Used internally in the `MixAll' package.
#'
#' @keywords internal
#' @rdname initialize-methods
setMethod(
    f="initialize",
    signature=c("LearnAlgo"),
    definition=function(.Object,algo,nbIteration,epsilon)
    {
      # for algo
      if(missing(algo)) {.Object@algo<-"SIMUL"}
      else  {.Object@algo<-algo}
      # for epsilon
      if( missing(epsilon) ){ .Object@epsilon<-1e-07 }
      else{.Object@epsilon<-epsilon}
      # for nbIteration
      if(missing(nbIteration)){ .Object@nbIteration<-200 }
      else{.Object@nbIteration<-nbIteration}
      # validate
      validObject(.Object)
      return(.Object)
    }
)

#' @rdname print-methods
#' @aliases print print-algo,LearnAlgo-method
setMethod(
    f="print",
    signature=c("LearnAlgo"),
    function(x,...){
      cat("****************************************\n")
      cat("*** MixAll LearnAlgo:\n")
      cat("* algorithm            = ", x@algo, "\n")
      cat("* number of iterations = ", x@nbIteration, "\n")
      cat("* epsilon              = ", x@epsilon, "\n")
      cat("****************************************\n")
    }
)

#' @rdname show-methods
#' @aliases show show-algo,LearnAlgo-method
setMethod(
    f="show",
    signature=c("LearnAlgo"),
    function(object){
      cat("****************************************\n")
      cat("*** MixAll LearnAlgo:\n")
      cat("* algorithm            = ", object@algo, "\n")
      cat("* number of iterations = ", object@nbIteration, "\n")
      cat("* epsilon              = ", object@epsilon, "\n")
      cat("****************************************\n")
    }
)

#' @rdname extract-methods
#' @aliases [,LearnAlgo-method
setMethod(
  f="[",
  signature(x = "LearnAlgo"),
  definition=function(x, i, j, drop){
    if ( missing(j) ){
      switch(EXPR=i,
          "algo"={return(x@algo)},
          "nbIteration"={return(x@nbIteration)},
          "epsilon"={return(x@epsilon)},
          stop("This attribute doesn't exist !")
      )
    }else{
      stop("This attribute is not a list !")
    }
  }
)

#' @name [
# @docType methods
#' @rdname extract-methods
#' @aliases [<-,LearnAlgo-method
setReplaceMethod(
    f="[",
    signature(x = "LearnAlgo"),
    definition=function(x,i,j,value){
      if ( missing(j) )
      {
        switch(EXPR=i,
            "algo"={x@algo<-value},
            "nbIteration"={x@nbIteration<-value},
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
