## #################################################################
##
## setEvidence etc.
## 
## The old setFinding type functions (in the Finding.R file) are just
## remapped.
##
## setEvidence with the syntax below is used in Scutari's bnlearn
## book.
##
## As of August 2016, these functions call a new implementation.
##
#################################################################

## FIXME: setEvidence: Could check that invalid nodes and/or invalid states are not given.
## FIXME: setEvidence: Functions .evidenceHasValidStates and .evidenceHasValidNodes
## FIXME: setEvidence: Does the checking. For now it has just been noted in the .Rd files
## FIXME: setEvidence: invalid states / nodes are ignored

#' @title Set, update and remove evidence.
#' 
#' @description Set, update and remove evidence.
#' 
#' @name grain_evidence
#' 
#' @aliases setEvidence retractEvidence absorbEvidence
#' @param object A "grain" object
#' @param nodes A vector of nodes; those nodes for which the
#'     (conditional) distribution is requested.
#' @param states A vector of states (of the nodes given by 'nodes')
#' @param evidence An alternative way of specifying findings
#'     (evidence), see examples below.
#' @param nslist deprecated
#' @param propagate Should the network be propagated?
#' @param details Debugging information
#'
#' @return A list of tables with potentials.
#'
#' @note \code{setEvidence()} is an improvement of \code{setFinding()}
#'     (and as such \code{setFinding} is obsolete). Users are
#'     recommended to use \code{setEvidence()} in the future.
#' 
#' \code{setEvidence()} allows to specification of "hard evidence" (specific
#' values for variables) and likelihood evidence (also known as virtual
#' evidence) for variables.
#' 
#' The syntax of \code{setEvidence()} may change in the future.
#'
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{setFinding}}, \code{\link{getFinding}},
#'     \code{\link{retractFinding}}, \code{\link{pFinding}}
#' @references Søren Højsgaard (2012). Graphical Independence Networks
#'     with the gRain Package for R. Journal of Statistical Software,
#'     46(10), 1-26.  \url{http://www.jstatsoft.org/v46/i10/}.
#' @keywords models utilities
#' @examples
#' 
#' data(chest_cpt)
#' chest.bn <- grain(compileCPT(chest_cpt))
#' chest.bn <- compile(chest.bn)
#' 
#' ## 1) These two forms are identical
#' setEvidence(chest.bn, c("asia", "xray"), c("yes", "yes"))
#' setFinding(chest.bn, c("asia", "xray"), c("yes", "yes"))
#' 
#' ## 2) Suppose we do not know with certainty whether a patient has
#' ## recently been to Asia. We can then introduce a new variable
#' ## "guess.asia" with "asia" as its only parent. Suppose
#' ## p(guess.asia=yes|asia=yes)=.8 and p(guess.asia=yes|asia=no)=.1
#' ## If the patient is e.g. unusually tanned we may set
#' ## guess.asia=yes and propagate.
#' ##
#' ## This corresponds to modifying the model by the likelihood (0.8,
#' ## 0.1) as
#' 
#' setEvidence(chest.bn, c("asia", "xray"), list(c(0.8, 0.1), "yes"))
#' 
#' ## 3) Hence, the same result as in 1) can be obtained with
#' setEvidence(chest.bn, c("asia", "xray"), list(c(1, 0), "yes"))
#' 
#' ## 4) An alternative specification using evidence is
#' setEvidence(chest.bn, evidence=list(asia=c(1, 0), xray="yes"))
#' 

#' @rdname grain_evidence
#' @export 
setEvidence <- function(object, nodes=NULL, states=NULL, evidence=NULL, nslist=NULL,
                        propagate=TRUE, details=0){
    
    if (!is.null(nslist))
        stop("Argument 'nslist' has been deprecated; please use 'evidence' instead\n")

    if ( is.null( evidence ) && is.null( nodes ) )
        stop( "Evidence is not given; nothing to do...")
    
    if ( is.null( evidence ) ) ## Then 'nodes' must be given
        evidence <- .nodes.states2evidence( nodes, states )      
    
    setEvi_(object, evidence, propagate=propagate, details=details)
    
}

#' @rdname grain_evidence
#' @export 
retractEvidence <- function(object, nodes=NULL, propagate=TRUE){
    retractEvi_(object, items=nodes, propagate=propagate)
}

#' @export 
#' @rdname grain_evidence
absorbEvidence <- function(object, propagate=TRUE ){
    absorbEvi( object, propagate=propagate )
    ## FIXME Should be absorbEvi_ 
}


setEvi <- function(object, nodes=NULL, states=NULL, evidence=NULL, 
                   propagate=TRUE, details=0){
    if ( !inherits(object, "grain") )
        stop("'object' is not a 'grain' object")
    
    if ( is.null( evidence ) && is.null( nodes ) )
        stop( "Evidence is not given; nothing to do...")
    
    if ( is.null( evidence ) ) ## Then 'nodes' must be given
        evidence <- .nodes.states2evidence( nodes, states )      
    
    setEvi_( object, evidence=evidence, propagate=propagate, details=details)
}

setEvi_ <- function(object, evidence=NULL, propagate=TRUE, details=0){
#    details=0
#    cat("++++ setEvi_ input evidence: \n"); str(evidence)

    ## If object is not compiled then do so
    
    if (is.null_ev( evidence )){
        cat("Nothing to do\n")
    } else {
        if (!isCompiled(object)){
            object <- compile(object)
        } 
        
        oe <- getEvidence( object ) # Der er noget med typen (old/new)
        if (details > 0){ print("old.evidence"); print(oe) }
        ne <- new_ev( evidence, universe(object)$levels ) 
        if (details > 0){ print("new evidence"); print( ne ) }
        
        ## Hvis der er eksisterende evidens, så skal det tages ud af det nye
        if (!is.null_ev( oe ) ){
            ne <- setdiff_ev( ne, oe )
            if (details>0) { print("new evidence - after modif"); print( ne ) }
        }

        if (length(varNames(ne)) > 0){
            rp  <- getgrain(object, "rip")    
            host  <- getHostClique(varNames(ne), rp$cliques)
            #str(list(ne, host))
            object$potential$pot_temp <- insertEvi(ne, getgrain(object, "pot_temp"), host)
            
            te <- if (is.null_ev(oe)) ne else union_ev(oe, ne)
            ##te <- union_ev( oe, ne )
            if (details>0) {print("total evidence"); print( te )}
            object$evidence <- te
        }         
    } 
    if (propagate) propagate(object) else object
}


retractEvi <- function(object, items=NULL, propagate=TRUE){
    if ( !inherits(object, "grain") )
        stop("'object' is not a 'grain' object")
    retractEvi_(object, items=items, propagate=propagate)
}

## #' @rdname grain-evi
retractEvi_ <- function(object, items=NULL, propagate=TRUE){
    ##cat("++++ retractEvidence_\n")
    .resetgrain <- function(x){
        ##x$temppot       <- x$pot_orig
        ##x$potential$pot_temp       <- pot(x)$pot_orig
        x$potential$pot_temp       <- getgrain(x, "pot_orig")
        x$evidence       <- NULL
        isPropagated(x) <- FALSE
        x
    }

    if ( is.null( items ) ){
        object <- .resetgrain( object )
    } else {
        if (!(is.character( items ) || is.numeric( items )) )
            stop("'items' must be a character or numeric vector")
        oe <- getEvidence( object )
        if (!is.null_ev( oe )){
            vn <- varNames( oe )
            if (is.numeric(items)){
                items <- vn[items]
            }
            keep <- setdiff( vn, items)
            ## NB: keep <= varNames            
            ## hvis keep==varNames(oe) så gør intet
            ## hvis keep=Ø så bare reset
            ## hvis Ø < keep < varNames så gør som nedenfor
            if ( length( keep ) < length( vn ) ){
                object <- .resetgrain( object )                
                if (length( keep ) > 0){
                    ne <- subset( oe, select=keep )
                    object <- setEvi( object, evidence = ne )
                }
            }
        }
    }
    if (propagate) propagate(object) else object
}

## #' @name grain-evi
absorbEvi <- function(object, propagate=TRUE ){
    if ( !inherits(object, "grain") )
        stop("'object' is not a 'grain' object")
    absorbEvi_( object, propagate = propagate )
}

## #' @name grain-evi
absorbEvi_<- function(object, propagate=TRUE ){
    ## Update 'object' as
    ## 1) set pot_orig <- pot_temp
    ## 2) ignore any finding
    object$potential$pot_orig <-  ## FRAGILE assignment
        getgrain(object, "pot_temp")
    object$evidence <- NULL
    isPropagated(object) <- FALSE

    if (propagate) propagate(object) else object
}

#' @export 
#' @name grain_evidence
pEvidence <- function(object){
    if ( !inherits(object, "grain") )
        stop("'object' is not a 'grain' object")
                                        #attr(pot(object)$pot_equi, "pEvidence")
    attr(getgrain(object, "pot_equi"), "pEvidence")
}

## #' @name grain-evi
pEvi <- function(object)
    pEvidence(object)

#' @export 
#' @name grain_evidence
getEvidence <- function(object){
    if ( !inherits(object, "grain") )
        stop("'object' is not a 'grain' object")
    object$evidence
}

## #' @name grain-evi
"dropEvi<-" <- function(object, value=NULL){
    retractEvi(object, value)
}

## #' @name grain-evi
"addEvi<-" <- function(object, value=NULL){
    setEvi(object, value)
}

#' @export 
#' @name grain_evidence
evidence <- function(object){
    UseMethod("evidence")
}

#' @export 
#' @name grain_evidence
evidence.grain <- function(object){
    getEvidence( object )
}

#' @name grain_evidence
#' @param value The evidence in the form of a named list or an evidence-object. 
"evidence<-" <- function(object, value=NULL){
    UseMethod("evidence<-")
}

#' @name grain_evidence
"evidence<-.grain" <- function(object, value=NULL){
    if (is.null( value ))
        retractEvi( object )
    else
        setEvidence(object, evidence=value)
}

## Alternative names

## #' @name grain-evi
addEvi  <- setEvi

## #' @name grain-evi
dropEvi <- retractEvi

## #' @name grain-evi
getEvi  <- getEvidence

## #############################################################
##
## UTILITIES
##
## #############################################################

insertEvi <- function(evi.list, pot, hostclique){
    ##cat("insertEvi\n")
    if ( !inherits(evi.list, "grain_ev") )
        stop("'object' is not a 'grain_ev' object")
    
    for (i in seq_along( evi.list$evidence ) ){
        p <- evi.list$evidence[[ i ]]
        #print(p)
        j <- hostclique[ i ]
        #print( j )
        pot[[j]] <- tabMult( pot[[ j ]], p )
    }
    pot
}


.nodes.states2evidence <- function(nodes, states){
    if (!is.null( states ) && length( nodes )==length( states )){
        evidence <- as.vector(states, "list")
        names(evidence) <- nodes
        evidence
    } else {
        stop( "Either states are not given or nodes and states do not have same length" )
    }
}

getHostClique <- function(set.list, cliques){
    out <- lapply(set.list,
                  function(x){
                      ##get_superset_(x, cliques, all=FALSE)
                      get_superset(x, cliques, all=FALSE)
                  })
    len <- sapply(out, length)
    if (any( len == 0)){
        message("Set(s) not in any clique:")
        str(set.list[ len==0 ])
        stop("exiting...\n")
    }
    unlist( out )
}






