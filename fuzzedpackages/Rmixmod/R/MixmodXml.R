##################################################################################
#                               MixmodXml.R                                     ##
##################################################################################

#' @include MixmodCluster.R
#' @include MixmodLearn.R
#' @include MixmodPredict.R
NULL

#' Constructor of [\code{\linkS4class{MixmodXmlCheck}}] class
#' 
#' This is a class to handle XML files (TODO: describe...)
#'
#' @name MixmodXmlCheck-class
#' @rdname MixmodXmlCheck-class
#' @exportClass MixmodXmlCheck
#'
setClass(
  Class="MixmodXmlCheck",
  representation=representation(
    xmlFile = "character",
    xmlType = "character"
    )
  )
setMethod(
  f="initialize",
  signature=c("MixmodXmlCheck"),
  definition=function(.Object, xmlFile){
    .Object@xmlFile <- xmlFile
    .Object@xmlType <- "unknown"    
    return(.Object)
    }
  )

#' mixmodXmlCheck
#'
#' TODO: describe...
#'
#' @param ... ...
#'
#' @return Object of type MixmodXmlCheck
#'
#' @export
mixmodXmlCheck <- function(...){
  return(new("MixmodXmlCheck", ...))
}

#' mixmodXmlLoad
#'
#' TODO: describe...
#'
#' @param xmlFile ...
#' @param numFormat ...
#'
#' @return XML output of mixmod methods
#'
#' @export
mixmodXmlLoad <- function(xmlFile, numFormat="humanReadable"){
  xmlIn <- mixmodXmlInput(xmlFile, numFormat=numFormat, conversionOnly=TRUE)
  xem <- new("MixmodXmlCheck", xmlFile)
  .Call("xMain", xem, PACKAGE="Rmixmod")
  if(xem@xmlType == "clustering"){
    return(mixmodCluster(xmlIn=xmlIn))
  }
  if(xem@xmlType == "learn"){
    return(mixmodLearn(xmlIn=xmlIn))
  }
  if(xem@xmlType == "predict"){
    return(mixmodPredict(xmlIn=xmlIn))
  }
}
