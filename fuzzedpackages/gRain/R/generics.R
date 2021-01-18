#' @title gRain generics
#'
#' @description Generic functions etc for the gRain package
#'
#' @name generics
#' 
#' @param x,object A relevant object.
#' @param nodes Some nodes of the object.
#' @param ... Additional arguments; currently not used.
#' 

#' @rdname generics
nodeNames  <- function(x) UseMethod("nodeNames")

#' @rdname generics
nodeNames.grain  <- function(x)
  getgrain(x, "universe")$nodes

#' @rdname generics
nodeStates <- function(x, nodes=nodeNames(x)) UseMethod("nodeStates")

#' @rdname generics
nodeStates.grain <- function(x, nodes=nodeNames(x)){
  getgrain(x, "universe")$levels[nodes]
}

#' @rdname generics
universe <- function(object, ...) UseMethod("universe")

#' @rdname generics
universe.grain <- function(object, ...)
    getgrain(object, "universe")

#' @rdname generics
varNames.grainEvidence_ <- function(x)
    getgrain(x, "summary")$nodes

#' @rdname generics
rip.grain <- function(object, ...)
    getgin(object, "rip")



#' @rdname generics
vpar.cpt_spec <- function(object, ...){
    lapply(object, function(u) names(dimnames(u)))
}

#' @rdname generics 
vpar.cpt_grain <- function(object, ...){
    lapply(getgin(object, "cptlist"), function(u) names(dimnames(u)))
}

isCompiled <- function(x) getgin(x, "isCompiled")

isPropagated <- function(x) getgin(x, "isPropagated")

"isCompiled<-" <- function(object, value){
    object$isCompiled <- value
    object
}

"isPropagated<-" <- function(object, value){
    object$isPropagated <- value
    object
}




## #' @rdname generics
## uni <- function(object)
##     UseMethod("uni")

## #' @rdname generics
## uni.grain <- function(object)
##     getgin(object, "universe")

## #' @rdname generics
## pot <- function(object)
##     UseMethod("pot")

## #' @rdname generics
## pot.grain <- function(object)
##     getgin(object, "potential")


## #' @rdname generics
## cpt <- function(object)
##     UseMethod("cpt")

## #' @rdname generics
## cpt.cpt_grain <- function(object)
##     getgin(object, "cptlist")

## #' @rdname generics
## #' @param position Where to insert 'value'
## #' @param value Value to insert at 'position'
## "cpt<-" <- function(object, position, value){
##     UseMethod("cpt<-")
## }

## #' @rdname generics
## "cpt<-.cpt_grain" <- function(object, position, value){
##     object$cptlist[position] <- value
## }


## #' @rdname generics
## potential <- function(object)
##    UseMethod("potential")

## #' @rdname generics
## potential.grain <- function(object)
##    object$potential

