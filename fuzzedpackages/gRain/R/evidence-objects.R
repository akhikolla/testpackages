## ###############################################################
##
#' @title Evidence objects
#' @description Functions for defining and manipulating evidence.
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @name evidence_object
##
## ###############################################################
#'
#' @aliases subset.grain_ev print.grain_ev varNames.grain_ev
#' 
#' @details Evidence is specified as a list. Internally, evidence is
#'     represented as a grain evidence object which is a list with 4 elements.
#' 
#' @examples
#'
#' ## Define the universe
#' 
#' uni <- list(asia = c("yes", "no"), tub = c("yes", "no"), smoke = c("yes", "no"),
#'             lung = c("yes", "no"), bronc = c("yes", "no"), either = c("yes", "no"),
#'             xray = c("yes", "no"), dysp = c("yes", "no"))
#'
#' e1 <- list(dysp="no", xray="no")
#' eo1 <- new_ev( e1, levels=uni )
#' eo1
#' as.data.frame( eo1 )
#' eo1 %>% str
#' 
#' e1.2 <- list(dysp="no", xray=c(0, 1))
#' eo1.2 <- new_ev( e1.2, levels=uni )
#' eo1.2
#'
#' # Notice that in eo1.2, xray is not regarded as hard
#' # evidence but as a weight on each level. Other than that, eo1.2
#' # and eo1 are equivalent here. This is used in connection
#' # with specifying likelihood evidence. 
#' 
#' e2 <- list(dysp="yes", asia="yes")
#' eo2 <- new_ev(e2, uni)
#'
#' # If evidence 'e1' is already set in the network and new evidence
#' # 'e2' emerges, the evidence in the network must be updated. But
#' # there is a conflict in that dysp="yes" in 'e1' and
#' # dysp="no" in 'e2'. The (arbitrary) convention is that
#' # existsting evidence overrides new evidence so that the only new
#' # evidence in 'e2' is really asia="yes".
#'
#' # To subtract existing evidence from new evidence we can do:
#' setdiff_ev( eo2, eo1 )
#'
#' # Likewise the 'union' is
#' union_ev( eo2, eo1 )
#'
#' @export 
#' @rdname evidence_object
#' @param evi.list A named list with evidence; see 'examples' below.
#' @param levels A named list with the levels of all variables. 
new_ev <- function(evi.list=NULL, levels){

    if (inherits( evi.list, "grain_ev") ) return( evi.list )
    
    out <- list(nodes=character(0),
                is.hard.evidence=logical(0),
                hard.state=character(0),
                evidence=list() )
    
    if (length(evi.list) != 0){
        ## First remove all evidence specified as NA
        not.na <- !unlist(lapply(lapply(evi.list,is.na), any), use.names=FALSE)
        if (length( not.na ) > 0)
            evi.list <- evi.list[ not.na ]
        
        evidence           <- vector("list", length(evi.list))
        is.hard.evidence   <- rep.int(TRUE,  length(evi.list))
        hard.state         <- rep.int(NA,    length(evi.list))
        
        ## cat("new_ev\n"); print(evi.list)
        
        for (i in seq_along(evi.list)){
            ev <- evi.list[i]
            v <- ev[[1]]
            
            if( is.array(v)){
                n <- names(dimnames(v))
                is.hard.evidence[i]    <- FALSE
                evidence[[i]] <- v
                next
            }
            
            if (is.character(v)){
                n <- names(evi.list)[i]
                hard.state[i]  <- v
                evidence[[i]]  <- .hard.state2parray(n, v, levels[[n]])
                next
            }
            
            if (is.numeric(v)){
                n <- names(evi.list)[i]
                is.hard.evidence[i] <- FALSE
                evidence[[i]] <- .soft.state2parray(n, v, levels[[n]])
            }
        }
        
        ## If evidence is zero on all states or negative on some (or all) states then it is invalid
        keep <- unlist(lapply(evidence, function(e){ sum(e) !=0 && all(e>=0) }), use.names=FALSE)
        ## print(keep)
        
        nodes <- unique.default( unlist(lapply(evidence, .namesDimnames)),
                                use.names=FALSE )
        
        out$nodes            = nodes[keep]
        out$is.hard.evidence = is.hard.evidence[keep]
        out$hard.state       = hard.state[keep]
        out$evidence         = evidence[keep]
    }
    class( out ) <- "grain_ev"
    out
}

#' @export 
#' @rdname evidence_object
#' @param object Some R object.
is.null_ev <- function(object){
    if (missing(object)) TRUE
    else if (length(object)==0) TRUE
    else if (class(object)=="grain_ev" && length(varNames(object))==0) TRUE
    else FALSE

}

## #' @rdname evidence_object
## #' @param x Evidence object

#' @export
print.grain_ev <- function(x, ...){
    print( as.data.frame(x[1:3]) )
}


## #' @rdname evidence_object

varNames.grain_ev <- function(x) x$nodes



#' @rdname evidence_object
#' @param row.names Not used.
#' @param optional Not used.
#' @param x An evidence object.
#' @param ... Not used.
#' @export 
as.data.frame.grain_ev <-
    function (x, row.names = NULL, optional = FALSE, ...) {
        is.atom <- sapply(x, is.atomic)
        atom <- x[is.atom]
        n.atom <- length(atom)
        out <- as.data.frame(atom)
        notatom <- x[!is.atom]
        n <- names(notatom)
        for (i in 1:length(notatom)){
            out[i+n.atom] <- notatom[i]
        }
        out
    }

#' @export 
#' @rdname evidence_object
#' @param ev1,ev2 Evidence.
setdiff_ev <- function(ev1, ev2){
    if (length(ev1)==0) ev1 <- new_ev( ev1 )
    if (length(ev2)==0) ev2 <- new_ev( ev2 )
    nn  <- setdiff( varNames(ev1), varNames(ev2) )
    out <- subset(ev1, select=nn)
    class(out) <- "grain_ev"
    out
}

#' @export 
#' @rdname evidence_object
union_ev <- function(ev1, ev2){
    if (length(ev1)==0) ev1 <- new_ev( ev1 )
    if (length(ev2)==0) ev2 <- new_ev( ev2 )
    ev <- setdiff_ev( ev1, ev2 )
    out <- mapply(function(l1, l2){c(l1,l2)},
                  ev, ev2, SIMPLIFY=FALSE, USE.NAMES=TRUE)
    class(out) <- "grain_ev"
    out
}

#' @export 
subset.grain_ev <- function(x, subset, select, drop = FALSE, ...){
    if (missing(select)) x
    else if (length(select)==0) new_ev(list())
    else {
        nl <- as.list(1L:length(varNames(x)))
        names(nl) <- varNames(x)
        vars <- eval(substitute(select), nl, parent.frame())
        if (is.character(vars))
            vars <- match(vars, varNames(x))
        if (any(is.na(vars)))
            stop("'vars' contain NA")
        if (max(vars) > length(varNames(x)))
            stop("'vars' too large")
        out <- lapply(x, "[", vars)
        class(out) <- class(x)
        out
    }
}


## ###############################################
##
## UTILITIES
##
## ###############################################

## Bruges i new_ev
.hard.state2parray <- function(n, v, lev){
    #str(list(n,v,lev))
    tab <- .fast.parray(n, list(lev), rep.int(0, length(lev)))
    tab[ match( v, lev )] <- 1
    tab
}

## Bruges i new_ev
.soft.state2parray <- function(n, v, lev){
    #str(list(n,v,lev))
    .fast.parray(n, list(lev), v)
}

## Bruges af ovenstående fns
.fast.parray <- function(varNames, levels, values=1){
    #str(list(varNames, levels, values))
    dn <- if(is.list(levels)) levels else list(levels)
    names(dn) <- varNames
    array(values, dimnames=dn)
}
