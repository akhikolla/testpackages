
randomCPT <- function(object, states=c("yes", "no")){

    if (!inherits(object, "graphNEL"))
        stop("'object' must be a graphNEL\n")
    if (!is_dag(object))
        stop("'object' is not a DAG\n")

    vpa <- vpar( object )
    n.states  <- length(states)
    cpt <- lapply(vpa, function(zz)
                  cptable(zz, values=runif( n.states^length(zz) ), levels=states))

    compileCPT( cpt )
}

.print_probability <- function(vn){
    if (length(vn) > 1){
        cat(paste(" P(", vn[1], "|", paste(vn[-1], collapse=' '), ")\n"))
    } else {
        cat(paste(" P(", vn, ")\n"))
    }    
}

.formula2char <- function(f) {
    unlist(rhsf2list(f))
}

.namesDimnames <- function(x)
    names(dimnames(x))

setSliceValue <- function(x, slice, complement=FALSE, value=0){
    margin <- names(slice)
    level  <- unlist(slice, use.names=FALSE)
    idx <- tableGetSliceIndex(x, margin = margin, level = level,
                              complement = complement)
    x[idx] <- value
    x
}


listify_dots <- function(args){
    args <- lapply(args, function(a) if (!is.list(a)) list(a) else a)
    unlist(args, recursive=FALSE)    
}


#' @export 
getgrain<- function(object, name=c("universe", "data", "dag", "ug", "cptlist",
                                   "origpot", "temppot", "equipot",
                                   "pot_orig", "pot_temp", "pot_equi",
                                   "rip",
                                   "isCompiled", "isPropagated",
                                   "evidence", "pEvidence",
                                   "control", "details")){

    switch(name,
           universe 	    = object$universe,
           data 	    = object$data,
           dag 		    = object$dag,
           ug 		    = object$ug,
           cptlist          = object$cptlist,
           cpt              = object$cptlist,

           potential        = object$potential,
           origpot          = object$potential$pot_orig,
           temppot          = object$potential$pot_temp,
           equipot          = object$potential$pot_equi,

           pot_orig         = object$potential$pot_orig,
           pot_temp         = object$potential$pot_temp,
           pot_equi         = object$potential$pot_equi,

           rip              = object$rip,

           isCompiled       = object$isCompiled,
           isPropagated     = object$isPropagated,
           
           evidence         = object$evidence,
           pEvidence        = object$pEvidence,
           control          = object$control,
           details          = object$details
           )
}

#' @export 
getgin <- getgrain







.infoPrint <- function(details, limit=1, ...,  prefix='.'){
  if(details>=limit){
    cat(paste(rep(prefix, limit), collapse=''), "")
    cat(...)
  }
}

.infoPrint <- function(details, limit=1, ...,  prefix='.'){}


.infoPrint2 <- function(details, limit=1, fmt, ...,  prefix='.'){
  if (details>=limit)
    cat(paste(paste(rep(prefix, limit), collapse=''),  sprintf(fmt, ...), collapse=' '))
}

.colstr <- function(x, collapse=" ")
  paste(x, collapse=collapse)



.printList <- function(x){
  mapply(function(xx,ii) cat(" ",ii,paste(xx, collapse=' '),"\n"), x, 1:length(x))
  return()
}


#' @export
printlist <- function(x,d=0) UseMethod("printlist")

#' @export
printlist.default <- function(x,d=0){
  paste("(", paste(x,collapse=' '),")",sep='')
}

#' @export
printlist.list <- function(x,d=0){
  tmp     <- unlist(lapply(x, printlist, d+2),recursive=FALSE)
  prefix  <- as.list(c("(",rep(" ",length(tmp)-1)))
  posfix  <- as.list(c(rep(" ",length(tmp)-1),")"))
  as.list(mapply(function(l,x,r) {paste(l,x,r,sep='')}, prefix, tmp, posfix))
}

#' @export
splitVec <- function(val, lev) UseMethod("splitVec")

#' @export
splitVec.default <- function(val, lev){
  m    <- matrix(val,ncol=lev)
  cval <- unlist(apply(m,2,list),recursive=FALSE)
  cval
}

#' @export
splitVec.list <- function(val, lev){
  lapply(val, splitVec, lev)
}



## ###############################################################
##
## as.grain() methods
##
## added July 2014
##
## ###############################################################

## as.grain <- function(x, ...){
##     UseMethod("as.grain")
## }

## as.grain.cpt_spec <- function(x, ...){
##     grain( x )
## }

## as.grain.graphNEL <- function(x, data, smooth=0, ...){
##     if (missing(data))
##         stop("'data' needed to create network from graph")
##     grain(x, data=data, smooth=smooth)
## }
