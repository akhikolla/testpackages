##################################################################################
#                            MixmodPredict.R                                    ##
##################################################################################

#' @include global.R
#' @include MixmodResults.R
#' @include Mixmod.R
NULL

#' Constructor of [\code{\linkS4class{MixmodPredict}}] class
#'
#' This is a class to run discriminant analysis with mixmod.
#' 
#' \describe{
#'   \item{data}{numeric vector, matrix, or data frame of observations. Either qualitative or quantitative.}
#'   \item{dataType}{character. It defines whether data are quantitative or qualitative.}
#'   \item{nbVariable}{integer. The number of variables.}
#'   \item{nbSample}{integer. The number of observations.}
#'   \item{error}{a character. The mixmod error.}
#'   \item{classificationRule}{a [\code{\linkS4class{MixmodResults}}] object containing the classification rule.}
#'   \item{partition}{a matrix containing observations to predict.}
#'   \item{proba}{a matrix of probabilities.}
#' }
#'
#' @examples
#'   # start by extract 10 observations from iris data set
#'   remaining.obs<-sample(1:nrow(iris),10)
#'
#'   # then run a mixmodLearn() analysis without those 10 observations
#'   learn<-mixmodLearn(iris[-remaining.obs,1:4], iris$Species[-remaining.obs])
#'   # create a MixmodPredict to predict those 10 observations
#'   new("MixmodPredict", data=iris[remaining.obs,1:4], classificationRule=learn["bestResult"])
#'
#'   getSlots("MixmodPredict")
#'
#' @name MixmodPredict-class
#' @rdname MixmodPredict-class
#' @exportClass MixmodPredict
#'
setClass(
  Class="MixmodPredict",
  representation=representation(
    data = "matrix",
    dataType = "character",
    nbVariable = "integer",
    nbSample = "integer",
    error = "character",
    classificationRule = "MixmodResults",
    partition = "integer",
    proba = "matrix",
    xmlIn = "MixmodXmlInput",
    xmlOut = "character",
    trace = "numeric",
    massiccc = "numeric"        
  ),
  prototype=prototype(
    data = matrix(nrow=0,ncol=0),
    dataType = character(0),
    nbVariable = integer(0),
    nbSample = integer(0),
    error = "No error",
    partition = integer(0),
    proba = matrix(nrow=0,ncol=0),
    xmlOut = character(0),
    trace = numeric(0),
    massiccc = numeric(0)        
  ),
  # define validity function
  validity=function(object){
    # check for missing values
    if ( sum(is.na(object@data)) ){
      stop("data contains NA. Mixmod cannot deal with missing values!")
    }
    return(TRUE)
  }
)

#' Create an instance of the [\code{\linkS4class{MixmodPredict}}] class
#'
#' This function computes the second step of a discriminant analysis. The aim of this step is to assign remaining observations to one of the groups.
#' 
#' @param data matrix or data frame containing quantitative,qualitative or composite data. Rows correspond to observations and columns correspond to variables.
#' @param classificationRule a [\code{\linkS4class{MixmodResults}}] object which contains the classification rule computed in the mixmodLearn() or mixmodCluster() step.
#' @param ... ...
#'
#' @examples
#'
#'   # start by extract 10 observations from iris data set
#'   remaining.obs<-sample(1:nrow(iris),10)
#'   # then run a mixmodLearn() analysis without those 10 observations
#'   learn<-mixmodLearn(iris[-remaining.obs,1:4], iris$Species[-remaining.obs])
#'   # create a MixmodPredict to predict those 10 observations
#'   prediction <- mixmodPredict(data=iris[remaining.obs,1:4], classificationRule=learn["bestResult"])
#'   # show results
#'   prediction
#'   # compare prediction with real results
#'   paste("accuracy= ",mean(as.integer(iris$Species[remaining.obs]) == prediction["partition"])*100
#'         ,"%",sep="")
#' 
#'   ## A composite example with a heterogeneous data set
#'   data(heterodatatrain)
#'   ## Learning with training data
#'   learn <- mixmodLearn(heterodatatrain[-1],knownLabels=heterodatatrain$V1)
#'
#' @author Florent Langrognet and Remi Lebret and Christian Poli ans Serge Iovleff, with contributions from C. Biernacki and G. Celeux and G. Govaert \email{contact@@mixmod.org}
#' @return Returns an instance of the [\code{\linkS4class{MixmodPredict}}] class which contains predicted partition and probabilities.
#' @export
#'
oldmixmodPredict <- function(data, classificationRule, ...) {
  # non documented params
  dots <- list(...)
  dotnames<-as.list(names(dots))
  trace <- 0
  if("trace" %in% dotnames) trace <- dots["trace"]
  massiccc <- 0
  if("massiccc" %in% dotnames) massiccc <- dots["massiccc"]
  # check options
  
  # check whether there is classification rule
  if ( missing(classificationRule) ){
    stop("classificationRule is missing !")
  }

  if(missing(data)){
    stop("data is missing !")
  }
  if (!is.data.frame(data) & !is.vector(data) & !is.factor(data) ){
    stop("data must be a data.frame or a vector or a factor")
  }
  
  # create Mixmod object
  xem <- new( "MixmodPredict", data=data, classificationRule=classificationRule, trace=trace, massiccc=massiccc )
  # call predictMain
  .Call("predictMain", xem, PACKAGE="Rmixmod")
  # mixmod error?
  if ( xem@error != "No error" ) warning( paste("Mixmod error: ", xem@error) )
  
  # return MixmodPredict object
  return(xem)
}

#' mixmodPredict
#'
#' TODO: describe
#'
#' @param ... ...
#'
#' @return A MixmodPredict object
#'
mixmodPredict <- function(...) {
  # create Mixmod object
  xem <- new( "MixmodPredict", ...)
  # call predictMain
  .Call("predictMain", xem, PACKAGE="Rmixmod")
  # mixmod error?
  if ( xem@error != "No error" ) warning( paste("Mixmod error: ", xem@error) )
  return(xem)
}

#' Create an instance of the [\code{\linkS4class{MixmodPredict}}] class using new/initialize.
#' 
#' Initialization method. Used internally in the `Rmixmod' package.
#' 
#' @seealso \code{\link{initialize}}
#'
#' @keywords internal
#'
#' @rdname initialize-methods
#'
setMethod(
  f="initialize",
  signature=c("MixmodPredict"),
  definition=function(.Object,data,classificationRule,  xmlIn=NULL, xmlOut=NULL, trace=0, massiccc=0){
    if(!missing(xmlIn)){
        if(!missing(data)||!missing(classificationRule)){
          stop("xmlIn argument is mutually exclusive with all other arguments but xmlOut and trace");
        }
        .Object@xmlIn = xmlIn
        .Object@classificationRule <- new("MixmodResults")
        .Object@classificationRule@parameters = new("CompositeParameter")
      } else { # i.e without xmlIn
        
        if(!missing(data)){
          if (!is.data.frame(data) & !is.vector(data) & !is.factor(data) ){
            stop("data must be a data.frame or a vector or a factor")
          }
          if ( is.factor(data) ) data<-as.integer(data)
          else if ( !is.vector(data) ){
                                        # loop over columns to check whether type is factor
            for ( j in 1:ncol(data) ){
              if ( is.factor(data[,j]) ) data[,j] <- as.integer(data[,j])
            }
          }
          .Object@data <- as.matrix(data)
          .Object@nbSample <- nrow(.Object@data)
          .Object@nbVariable <- ncol(.Object@data)
        } else {
          stop("data is missing !")
        } 
                                        # check whether there is classification rule
        if ( missing(classificationRule) ){
          stop("classificationRule is missing !")
        }
        else if ( !is(classificationRule, "MixmodResults") ){
          stop("classificationRule must be a 'MixmodResults' object!")
        }
        else{
          .Object@classificationRule <- classificationRule
        }
        
        if ( is(classificationRule@parameters, "MultinomialParameter") ){
          .Object@dataType <- "qualitative"
        }
        else if ( is(classificationRule@parameters, "GaussianParameter") ){
          .Object@dataType <- "quantitative"
        } 
        else if ( is(classificationRule@parameters, "CompositeParameter") ){
          .Object@dataType <- "composite"
        }
        else{
          stop( "Unknown type of parameters within classificationRule!" )
        }
      } #end without xml
    .Object@trace = trace
    .Object@massiccc = massiccc    
    if(missing(xmlIn)){validObject(.Object)}
    return(.Object)
  }
)

#' @rdname print-methods
#' @aliases print print,MixmodPredict-method
#'
setMethod(
  f="print",
  signature=c("MixmodPredict"),
  function(x,...){
    cat("****************************************\n")
    cat("*** INPUT:\n")
    cat("****************************************\n")    
    print(x@classificationRule)
    if(length(x@data)!=0){
      cat("* data      =\n")
      print(x@data)
    }
    else{}
    if ( x@error == "No error" ){
      cat("\n\n")
      cat("****************************************\n")
      cat("*** PREDICTION:\n")
      cat("****************************************\n")
      cat("* partition     = ", x@partition, "\n" );
      if ( nrow(x@proba)>1 ){
        cat("* probabilities = |",formatC(x@proba[1,],digits=4,width=8,format="f"),"|\n")
        for ( i in 2:nrow(x@proba)){
          cat("                  |", formatC(x@proba[i,],digits=4,width=8,format="f"),"|\n")
        }
      }else{
        cat("* probabilities = ",formatC(x@proba,digits=4,format="f"),"\n")
      }
      cat("****************************************\n")
    }
    else{
      cat("\n\n")
      cat("****************************************\n")
      cat("*** NO OUTPUT - All models got errors !\n")
      cat("****************************************\n")
    }
  }
)

#' @rdname show-methods
#' @aliases show show,MixmodPredict-method
#'
setMethod(
  f="show",
  signature=c("MixmodPredict"),
  function(object){
    cat("****************************************\n")
    cat("*** INPUT:\n")
    cat("****************************************\n")
    print(object@classificationRule)
    if(length(object@data)!=0){
      nrowShow <- min(10,nrow(object@data))
      ncolShow <- min(10,ncol(object@data))
      cat("* data (limited to a 10x10 matrix) =\n")
      print(formatC(object@data[1:nrowShow,1:ncolShow]),quote=FALSE)
    }else{}
    cat("* ... ...\n")  
    
    if ( object@error == "No error" ){
      cat("\n\n")
      cat("****************************************\n")
      cat("*** PREDICTION:\n")
      cat("****************************************\n")
      cat("* partition     = ", object@partition, "\n" );
      if ( nrow(object@proba)>1 ){
        cat("* probabilities = |",formatC(object@proba[1,],digits=4,width=8,format="f"),"|\n")
        for ( i in 2:nrow(object@proba)){
          cat("                  |", formatC(object@proba[i,],digits=4,width=8,format="f"),"|\n")
        }
      }else{
        cat("* probabilities = ",formatC(object@proba,digits=4,format="f"),"\n")
      }
      cat("****************************************\n")
    }
    else{
      cat("\n\n")
      cat("****************************************\n")
      cat("*** NO OUTPUT - All models got errors !\n")
      cat("****************************************\n")
    }
  }
)

#' @rdname summary-methods
#' @aliases summary summary,MixmodPredict-method
#'
setMethod(
  f="summary",
  signature=c("MixmodPredict"),
  function(object, ...){
    if ( object@error == "No error" ){
      cat("**************************************************************\n")
      cat("* partition     = ", object@partition, "\n" );
      if ( nrow(object@proba)>1 ){
        cat("* probabilities = |",formatC(object@proba[1,],digits=4,width=8,format="f"),"|\n")
        for ( i in 2:nrow(object@proba)){
          cat("                  |", formatC(object@proba[i,],digits=4,width=8,format="f"),"|\n")
        }
      }else{
        cat("* probabilities = ",formatC(object@proba,digits=4,format="f"),"\n")
      }
      cat("**************************************************************\n")
    }
    return(invisible())
  }
)

#' @rdname extract-methods
#' @aliases [,MixmodPredict-method
#'
setMethod(
  f="[", 
  signature(x = "MixmodPredict"),
  definition=function (x, i, j, drop) {
    if ( missing(j) ){
      switch(EXPR=i,
        "data"={return(x@data)},
        "dataType"={return(x@dataType)},
        "nbCluster"={return(x@nbCluster)},
        "classificationRule"={return(x@classificationRule)},
        "partition"={return(x@partition)},
        "proba"={return(x@proba)},
        "error"={return(x@error)},
        stop("This attribute doesn't exist !")
      )
    }else{
      switch(EXPR=i,
        "data"={return(x@data[,j])},
        "partition"={return(x@partition[j])},
        "proba"={return(x@proba[,j])},
        stop("This attribute doesn't exist !")
      )
    }
  }
)
