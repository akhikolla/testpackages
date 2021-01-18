## FIXME: print.marg_spec missing
## FIXME: summary.pot_spec missing
## FIXME: summary.marg_spec missing

#' @title Compile conditional probability tables / cliques potentials.
#' 
#' @description Compile conditional probability tables / cliques
#'     potentials as a preprocessing step for creating a graphical
#'     independence network
#'
#' @name components_gather
#' 
#' @param x To \code{compileCPT} x is a list of conditional
#'     probability tables; to \code{compilePOT}, x is a list of clique
#'     potentials.
#'
## #' @param object A list of potentials or of CPTs.
#'
#' @param forceCheck Controls if consistency checks of the probability
#'     tables should be made.
#' 
#' @param ... Additional arguments; currently not used.
#' 
#' @aliases parse_cpt, parse_cpt.xtabs,parse_cpt.cptable_class, parse_cpt.default
#' 
#' @details
#'     * `compileCPT` is relevant for turning a collection of
#'     cptable's into an object from which a network can be built. For
#'     example, when specification of a cpt is made with cptable then
#'     the levels of the node is given but not the levels of the
#'     parents. `compileCPT` checks that the levels of variables in
#'     the cpt's are consistent and also that the specifications
#'     define a dag.
#' 
#'     * `compilePOT` is not of direct relevance for the
#'     user for the moment. However, the elements of the input should
#'     be arrays which define a chordal undirected graph and the
#'     arrays should, if multiplied, form a valid probability density.
#'  
#' @return A list with a class attribute.
#' 
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#'
#' @seealso \code{\link{extractCPT}}, \code{\link{extractPOT}}, \code{\link{extractMARG}}
#' 
#' @references Søren Højsgaard (2012). Graphical Independence Networks
#'     with the gRain Package for R. Journal of Statistical Software,
#'     46(10), 1-26.  \url{http://www.jstatsoft.org/v46/i10/}.
#'
#' @keywords utilities
#'
#' @examples
#'
#' data(chest_cpt)
#' x <- compileCPT(chest_cpt)
#' class(x)
#' grain(x)
#' 
#' ## FIXME: compileCPT/compilePOT examples missing.

## foo <- function(x, ..., z){
##     args <- c(list(x), list(...))
##     listify_dots(args)
## }
## listify_dots <- function(args){
##     args <- lapply(args, function(a) if (!is.list(a)) list(a) else a)
##     unlist(args, recursive=FALSE)    
## }
## #' @details \code{compileCPT}, \code{compilePOT} are wrappers for
## #'     \code{compileCPT} and \code{compilePOT} and are kept for
## #'     backward compatibility.
## #'


#' @rdname components_gather
#' @export
compileCPT <- function(x, ..., forceCheck=TRUE){
    args <- c(list(x), list(...))
    args <- listify_dots(args)
    .compileCPT(args, forceCheck=forceCheck)
}

.compileCPT <- function(x, forceCheck=TRUE){
    ## x: A list of cpts (arrays)

    type <- is.list(x) ##&& all(sapply(x, is.named.array))
    if (!type) stop("A list of named arrays is expected")        ## FIXME Can also be cptable_class...
    
    ## zz: Internal representation of cpts
    zz  <- lapply(x, parse_cpt)

    universe <- .create_universe(zz)

    ## Given node names; need to check that they are not replicated
    vn_given <- sapply(zz, "[[", "vnam")
    if (length(vn_given) != length(unique(vn_given)))
        stop("Some nodes specified more than once: ", toString(vn_given))
    
    ## Are all cpts defined?
    ss <- setdiff(unique(unlist(vn_given)),  universe$nodes)
    if (length(ss) > 0)
        stop(paste("Distribution not specified for nodes(s):", toString(ss)))
    
    ## Does specification define a DAG? If x is cpt_rep the answer is yes
    if (inherits(x, "cpt_rep")){
        graph <- attr(x, "graph")
    } else {
        vp <- lapply(zz, "[[", "vpar")
        graph <- dagList(vp, forceCheck=forceCheck)
    }
    ## Need list of cpts (each represented as an array)
    out <- lapply(seq_along(zz), .create_array, zz, universe)    
    names(out) <- universe$nodes

    attr(out, "universe") <- universe
    attr(out, "dag")    <- graph
    class(out)            <- "cpt_spec"
    out
}

#' @rdname components_gather
compilePOT <- function(x, ..., forceCheck=TRUE){
    args <- c(list(x), list(...))
    args <- listify_dots(args)
    .compilePOT(args, forceCheck=forceCheck)   
}


## #############################################################

#' @export
print.cpt_spec <- function(x, ...){
    cat("cpt_spec with probabilities:\n")
    lapply(x,
           function(xx){
               vn <- varNames(xx)
               .print_probability(vn)
           })
    invisible(x)
}

#' @export
print.pot_spec <- function(x, ...){
    cat("pot_spec with potentials:\n")
    lapply(x,
           function(xx){
               vn <- names(dimnames(xx))
               cat("(", paste(vn, collapse=' '),") \n")
           })    
    invisible(x)
}

summary.cpt_spec <- function(object, ...){
    cat("cpt_spec with probabilities:\n")
    lapply(object,
           function(xx){
               vn <- varNames(xx)
               .print_probability(vn)               
           })
    invisible(object)

}

## ###########################################################
## Helper functions  -- used only in grain-main.R
## ###########################################################
as_cpt_spec_simple <- function(x){
    z <- c(x)
    attr(z, "universe") <- attr(x, "universe")
    class(z) <- "cpt_spec_simple"
    z
}

print.cpt_spec_simple <- function(x,...){
    cat("cpt_spec_simple with probabilities:\n")
    lapply(x,
           function(xx){
               vn <- varNames(xx)
               .print_probability(vn)               
           })
  invisible(x)
}



## ###################################################
##
## dot functions below here
##
## ###################################################

.compilePOT <- function(x, ...){
    ## x: a list of arrays, and a rip attribute

    type <- is.list(x) && all(sapply(x, is.named.array))
    if (!type) stop("A list of named arrays is expected")    

    universe  <- .make.universe(x)
    
    if (inherits(x, "pot_rep")){ ## Result of extractPOT
        graph <- attr(x, "graph")
        rp    <- attr(x, "rip")
    } else {
        graph <- lapply(x, .namesDimnames)
        graph <- ug(graph)
        rp    <- rip(graph)
    }
    
    attr(x, "universe") <- universe
    attr(x, "ug")    <- graph 
    attr(x, "rip")   <- rp
    class(x) <- "pot_spec"
    x
}

.create_universe <- function(zz){
    vn <- unlist(lapply(zz, "[[", "vnam"))
    vl <- lapply(zz, "[[", "vlev")
    di <- unlist(lapply(vl, length))
    names(vl) <- vn        
    universe  <- list(nodes = vn, levels = vl, nlev = di)
    universe
}

.create_array <- function(i, zz, universe){
    cp <- zz[[i]]
    dn <- universe$levels[cp$vpar]
    di <- sapply(dn, length)
    val <- array(rep(1.0, prod(di)), dim=di, dimnames=dn)
    if (length(cp$values) > 0)
        val[] <- cp$values + cp$smooth
    val
}

.make.universe <- function(x){
    lll       <- unlist(lapply(x, dimnames), recursive=FALSE)
    nnn       <- names(lll)
    iii       <- match(unique(nnn), nnn)
    levels    <- lll[iii]
    vn        <- nnn[iii]
    di        <- c(lapply(levels, length), recursive=TRUE)
    names(di) <- vn
    universe  <- list(nodes = vn, levels = levels, nlev   = di)
    universe
}









## ##################################################################
##
## Extend compilation function 
##
## ##################################################################

## compile.cpt_rep <- function(object, ...)
##     compileCPT(object, ...)

## compile.pot_rep <- function(object, ...)
##     compilePOT(object, ...)

## #################################################################


## ##################################################################
##
## INTERNAL UTILITIES
##
## The .parse_cpt functions are used only in compileCPT
##
## ##################################################################

#' @rdname components_gather
#' @param xi cpt in some representation
#' @export
parse_cpt <- function(xi){
    UseMethod("parse_cpt")
}

#' @export
parse_cpt.xtabs <- function(xi){
    NextMethod("parse_cpt")
}

#' @export
parse_cpt.cptable_class <- function(xi){
    .parse_cpt_finalize(varNames(xi), valueLabels(xi)[[1]],
                        as.numeric(xi), attr(xi, "smooth"))
}

#' @export
parse_cpt.default <- function(xi){
    if (!is.named.array(xi)) stop("'xi' must be a named array")
    .parse_cpt_finalize(varNames(xi), valueLabels(xi)[[1]],
                        as.numeric(xi), 0)
}

.parse_cpt_finalize <- function(vpar, vlev, values, smooth){

    ## Normalization of CPTs happen here
    values <- matrix(values, nrow=length(vlev))
    s  <- colSums(values)
    for (j in 1:ncol(values)) values[, j] <- values[, j] / s[j]
    values <- as.numeric(values)

    out <- list(vnam=vpar[1], vlev=vlev, vpar=vpar, values=values,
                normalize="first", smooth=smooth)
    class(out) <- "cpt_generic"
    out    
}











## compilePOT <- function(x, ...){
##     ## x: a list of arrays, and a rip attribute
    
##     .make.universe <- function(x){
##         lll       <- unlist(lapply(x, dimnames), recursive=FALSE)
##         nnn       <- names(lll)
##         iii       <- match(unique(nnn), nnn)
##         levels    <- lll[iii]
##         vn        <- nnn[iii]
##         di        <- c(lapply(levels, length), recursive=TRUE)
##         names(di) <- vn
##         universe  <- list(nodes = vn, levels = levels, nlev   = di)
##         universe
##     }

##     if (!inherits(x, "pot_rep")) stop("can not compile 'x'\n")
##     if (is.null(attr(x, "rip"))) stop("no rip attribute; not a proper POT_spec object")

##     attr(x, "universe") <- .make.universe(x)
##     attr(x, "ug")       <- ug(attr(x, "rip")$cliques)
##     class(x) <- "pot_spec"
##     x
## }





