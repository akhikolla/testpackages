
make_cptlist <- function(object){
    if (!inherits(object, "grain"))
        stop("Object is not a 'grain' object\n")
    if (!isCompiled(object)) 
        object <- compile(object)
    dg <- ug2dag(getgin(object, "ug"))
    vp <- vpar(dg)
    lapply(vp, function(vv){
        z <- qgrain(object, nodes=vv, type="cond")
        tabPerm(z, vv)}
        )
}



## #' @title Recent additions to gRain
## #'
## #' @description Recent additions to gRain
## #' 
## #' @name recent
## #'
## #' @aliases rip<- rip<-.grain mkcptlist as.cptlist
## #'
## #' @param object An appropriate R object
## #' @param value A value to be assigned
## #' @param ... Additional arguments
## #' 

## "rip<-" <- function(object, value)
##     UseMethod("rip<-")

## "rip<-.grain" <- function(object, value){
##     if (!is.null(value) && !inherits(value, "ripOrder")) stop("Invalid 'value'\n")
##     object$rip <- value
##     object$isCompiled <- FALSE
##     object$isPropagated <- FALSE
##     object$potential <- NULL 
##     object
## }

## as.cptlist <- function(..., forceCheck=TRUE, details=0){
##     args <- list(...)
##     if (inherits(args, "list") && length(args) == 1)
##         compileCPT(args[[1]])
##     else
##         compileCPT(args)
## }

