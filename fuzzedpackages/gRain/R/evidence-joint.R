## ###############################################################
##
#' @title Set joint evidence in grain objects
#' @description Setting and removing joint evidence in grain objects. 
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' 
## ###############################################################
#'
#' @name grain_jevidence
#'
#' @param object A "grain" object.
#' @param evidence A list of evidence. Each element is a named array.
#' @param propagate Should evidence be absorbed once entered; defaults
#'     to TRUE.
#' @param details Amount of printing; for debugging.
#' @param ev A named list.
#' @param levels A named list.
#'
#' @note All the joint evidence functionality should be used with
#'     great care.
#'
#' @aliases print.grain_jev
#'
#' @examples
#'
#' data(chest_cpt)
#' chest.bn <- grain(compileCPT(chest_cpt))
#' chest.bn <- compile(chest.bn)
#' 
#' uni <- list(asia = c("yes", "no"), tub = c("yes", "no"),
#'             smoke = c("yes", "no"), lung = c("yes", "no"),
#'             bronc = c("yes", "no"), either = c("yes", "no"),
#'             xray = c("yes", "no"), dysp = c("yes", "no"))
#'
#' ev <- list(tabNew("asia", levels=uni, values=c(1,0)),
#'            tabNew("dysp", levels=uni, values=c(1,0)),
#'            tabNew(c("dysp","bronc"), levels=uni, values=c(.1, .2, .9, .8)) )
#'
#' chest.bn
#' chest.bn2 <- setJEvidence(chest.bn, evidence=ev)
#' chest.bn2
#'
#' # Notice: The evidence is defined on (subsets of) cliques of the junction tree
#' # and therefore evidence can readily be absorbed:
#' getgrain(chest.bn, "rip")$cliques  %>% str
#'
#' # On the other hand, below is evidence which is not defined cliques
#' # of the junction tree and therefore evidence can not easily be
#' # absorbed.  Hence this will fail:
#'
#' \dontrun{
#' ev.fail <- list(tab(c("dysp","smoke"), levels=uni, values=c(.1, .2, .9, .8)) )
#' setJEvidence(chest.bn, evidence=ev.fail)
#' }
#'
#' ## Evidence can be removed with
#'
#' retractJEvidence(chest.bn2)      ## All evidence removed.
#' retractJEvidence(chest.bn2, 0)   ## No evidence removed.
#' retractJEvidence(chest.bn2, 1:2) ## Evidence items 1 and 2 are removed.
#' 
#' # Setting additional joint evidence to an object where joint
#' # evidence already is set will cause an error. Hence this will fail:
#'
#' \dontrun{
#'   ev2 <- list(smoke="yes")
#'   setJEvidence(chest.bn2, evidence=ev2)
#' }
#'
#' ## Instead we can do
#' new.ev <- c(getEvidence(chest.bn2), list(smoke="yes"))
#' chest.bn
#' setJEvidence(chest.bn, evidence=new.ev)
#' 
#' ## Create joint evidence object:
#' yn <- c("yes", "no")
#' db <- parray(c("dysp", "bronc"), list(yn, yn), values=c(.1, .2, .9, .8))
#' db
#' ev   <- list(asia=c(1, 0), dysp="yes", db)
#'
#' jevi <- new_jev(ev, levels=uni)
#' jevi
#'
#' chest.bn3 <- setJEvidence(chest.bn, evidence=jevi)
#' evidence(chest.bn3)
#'




#' @export 
#' @rdname grain_jevidence
setJEvidence <- function(object, evidence=NULL, propagate=TRUE, details=0){
    
    setJEvidence_(object, evidence=evidence, propagate=propagate, details=details)
}


setJEvidence_<- function(object, evidence=NULL, propagate=TRUE, details=0){

    if (!inherits(object, "grain")) stop("'object' is not a grain object")
    
    if (!is.null( evidence )) {
        if ( length(getEvidence( object )) > 0 )
            stop("'object' already has evidence set; can not proceed\n")
        
        if ( !inherits( evidence, "grain_jev" ) )
            evidence <- new_jev( evidence, universe( object )$levels )
        
        if (!isCompiled(object))
            object <- compile( object )
        
        vn  <- sapply(evidence, varNames)    
        rp  <- getgrain(object, "rip")    
        hc  <- getHostClique(vn, rp$cliques)
        pot <- getgrain(object, "pot_temp")
        pot2 <- insertJEvidence(evidence, pot, hc)
        object$potential$pot_temp <- pot2
        object$evidence <- evidence
    }
    object <- if (propagate) propagate( object ) else object
    object
}

## #' @rdname grain_jevidence
## #' @param evi.list A "grain_jev" object.
## #' @param pot A list of clique potentials (a potential is an array).
## #' @param hostclique A numerical vector indicating in which element of
## #'     'pot' each eviendence item in 'evi.list' should be inserted in.

insertJEvidence <- function(evi.list, pot, hostclique){
    if ( !inherits(evi.list, "grain_jev") )
        stop("'object' is not a 'grain_jev' object")
    #ee <<- evi.list
    for (i in seq_along( evi.list ) ){
        p <- evi.list[[ i ]]
        j <- hostclique[ i ]
        pot[[j]] <- tabMult( pot[[ j ]], p )
    }
    pot
}



#' @export 
#' @rdname grain_jevidence
#' @param items Items in the evidence list to be removed. Here,
#'     \code{NULL} means remove everything, \code{0} means nothing is
#'     removed. Otherwise \code{items} is a numeric vector.
#' 
retractJEvidence <- function(object, items=NULL, propagate=TRUE, details=0){
    if (! (is.numeric(items) || is.null(items) ))
        stop("'items' must be  numeric or NULL")            

    ev <- getEvidence( object )

    if (length( ev ) > 0) {        
        if ( is.null( items ) )
            object$evidence <- NULL ## remove all
        else {
            if (length(items)==1 && items==0){
                # do nothing
            } else {                     
                items <- items[ items > 0]
                keep <- seq_along( ev )[ -items ]
                if (any(is.na(keep))) stop("'items' out of range")
                ev <- ev[ keep]
                
                object$evidence <-
                    if (length(ev) > 0) ev else NULL
            }
        }
    }
    object <- if (propagate) propagate( object ) else object
    object
}

#' @export 
#' @rdname grain_jevidence
new_jev <- function(ev, levels){
    if (inherits(ev, "grain_jev")) return( ev )

    if (!is.list(ev)) stop("'ev' must be a list")
    if (!is.list(levels)) stop("'levels' must be a list")

    evname <- names(ev)
    spec <- !unlist(lapply(ev, is.named.array), use.names=FALSE)
    ## when spec=TRUE, we need to convert to arrays
    spec <- which( spec )
    if ( length( spec ) > 0 ){
        new.out <-
            lapply( spec, function(i){
                e <- ev[[i]]
                if ( !(is.character(e) || is.numeric(e)) )
                    stop(" evidence must be a character or a numeric vector ")
                n   <- names(ev)[i]
                if ( is.na( match(n, names(levels)) ))
                    stop("name does not exist in 'levels'")
                lev <- levels[[n]]
                if (is.numeric( e )){
                    if ( length(e) != length(lev) )
                        stop("length of evidence item not valid")
                    out <- array(e, dim=length(lev), dimnames=levels[n])
                    out
                } else {
                    m <- match( e, lev )
                    if (any(is.na(m)))
                        stop("evidence item contains invalid value")
                    out <- array(rep.int(0, length(lev)), dim=length(lev), dimnames=levels[n])
                    out[match( e, lev )] <- 1
                    out
                }
                out
            }        
            )
        ev[ spec ] <- new.out
    }    
    names(ev) <- NULL
    class(ev) <- "grain_jev"
    ev
}

## #' @rdname grain_jevidence
## #' @param x A "grain_jev" object.
## #' @param ... Additional arguments; currently not used.
print.grain_jev <- function(x, ...){
    vn <- lapply(x, varNames)
    str(vn)
    ##str(x)
}

