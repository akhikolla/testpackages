#' @title Propagate a graphical independence network (a Bayesian network)
#' 
#' @description Propagation refers to calibrating the cliques of the
#'     junction tree so that the clique potentials are consistent on
#'     their intersections; refer to the reference below for details.
#'
#' @name grain_propagate
#' 
#' @details The \code{propagate} method invokes \code{propagateLS}
#'     which is a pure R implementation of the Lauritzen-Spiegelhalter
#'     algorithm. The c++ based version is several times faster than
#'     the purely R based version.
#'
#' 
#' @aliases propagate.grain propagateLS propagateLS__
#' @param object A grain object
#' @param details For debugging info
#' @param engine Either "R" or "cpp"; "cpp" is the default and the
#'     fastest.
#' @param ... Currently not used
#' @return A compiled and propagated grain object.
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{grain}}, \code{\link[gRbase]{compile}}
#' @references Søren Højsgaard (2012). Graphical Independence Networks
#'     with the gRain Package for R. Journal of Statistical Software,
#'     46(10), 1-26.  \url{http://www.jstatsoft.org/v46/i10/}.
#' @keywords utilities models
#' @examples
#' 
#' yn   <- c("yes","no")
#' a    <- cptable(~asia,              values=c(1,99), levels=yn)
#' t.a  <- cptable(~tub+asia,          values=c(5,95,1,99), levels=yn)
#' s    <- cptable(~smoke,             values=c(5,5), levels=yn)
#' l.s  <- cptable(~lung+smoke,        values=c(1,9,1,99), levels=yn)
#' b.s  <- cptable(~bronc+smoke,       values=c(6,4,3,7), levels=yn)
#' e.lt <- cptable(~either+lung+tub,   values=c(1,0,1,0,1,0,0,1), levels=yn)
#' x.e  <- cptable(~xray+either,       values=c(98,2,5,95), levels=yn)
#' d.be <- cptable(~dysp+bronc+either, values=c(9,1,7,3,8,2,1,9), levels=yn)
#' chest.cpt <- compileCPT(list(a, t.a, s, l.s, b.s, e.lt, x.e, d.be))
#' chest.bn  <- grain(chest.cpt)
#' bn1  <- compile(chest.bn, propagate=FALSE)
#' bn2  <- propagate(bn1)
#' bn3  <- compile(chest.bn, propagate=TRUE)
#'  
#' @export propagate.grain

## propagate.grain: equipot is updated after propagation on temppot
## such that equipot will contain the updated potentials.
## object$equipot <- propagateLS(object$temppot,
##                               rip=object$rip, initialize=TRUE, details=details)

#' @export 
propagate.grain <- function(object, details=object$details, engine="cpp", ...){

    t0 <- proc.time()
    engine <- match.arg(tolower(engine), c("r", "cpp"))
    
    propfun <- switch(engine,
                   "r"   = {propagateLS},
                   "cpp" = {propagateLS__})

    object$potential$pot_equi <- ## FRAGILE assignment
        propfun(getgrain(object, "pot_temp"), rip=rip(object))

    isPropagated(object) <- TRUE
    
    ## FIXME: propagate.grain : Looks strange
    if (!is.null((ev <- getEvidence(object)))){
        attr(ev, "pEvidence") <- pEvidence(object)
        object$evidence <- ev
    }
    
    .timing(" Time: propagation:", object$control, t0)
    object
}

## Lauritzen Spiegelhalter propagation
##

## Don't remember the idea behind the 'initialize' argument; should always be true


#' @rdname grain_propagate
#' @param cqpotList Clique potential list
#' @param rip A rip ordering
#' @param initialize Always true
#' @export
propagateLS <- function(cqpotList, rip, initialize=TRUE, details=0){
    #details=20
    #cat(".Propagating BN: [propagateLS]\n")
    .infoPrint(details, 1, cat(".Propagating BN: [propagateLS]\n"))

    cliq       <- rip$cliques
    seps       <- rip$separators
    pa         <- rip$parent
    childList  <- rip$childList

    ncliq      <- length(cliq)

    ## This assignment is needed because RIP now returns 0 as the
    ## parent index for the first clique
    pa[pa==0]<-NA

    ## Backward propagation (collect evidence) towards root of junction tree
    ##
    .infoPrint(details,2, cat("..BACKWARD:\n"))
    t0 <- proc.time()
    if (ncliq > 1){
        for (i in ncliq:2){
            cq   <- cliq[[ i ]]
            sp   <- seps[[ i ]]
            .infoPrint2(details, 2, "Clique %d: {%s}\n",  i , .colstr( cq ))
            cq.pot   <- cqpotList[[ i ]]
            pa.pot   <- cqpotList[[pa[ i ]]]
            
            if (length(sp) >= 1 && !is.na(sp)){
                .infoPrint2(details, 2, "Marg onto sep {%s}\n", .colstr(sp))
                sp.pot               <- tableMargin(cq.pot, sp)
                cqpotList[[ i ]]     <- tableOp2(cq.pot, sp.pot, `/`)
                cqpotList[[pa[ i ]]] <- tableOp2(pa.pot, sp.pot, `*`)
            } else{
                zzz               <- sum(cq.pot)
                cqpotList[[1]]    <- cqpotList[[1]] * zzz
                cqpotList[[ i ]]  <- cq.pot / zzz
            }
        }
    }

    tmpd <- cqpotList[[1]]
    ##print(tmpd)
    normConst      <- sum(tmpd)
    tmpd <- tmpd / normConst
    cqpotList[[1]] <- tmpd
    ##print(tmpd)
    

    
    ## Forward propagation (distribute evidence) away from root of junction tree
    ##
    .infoPrint(details, 2, cat("..FORWARD:\n"))
    t0 <- proc.time()
    for ( i  in 1:ncliq){
        .infoPrint2(details, 2, "Clique %d: {%s}\n",  i , .colstr(cliq[[ i ]]))
        ch <- childList[[ i ]]
        if (length(ch) > 0)
        {
            .infoPrint2(details,2, "..Children: %s\n", .colstr(ch))
            for ( j  in 1:length(ch))
            {
                if (length(seps[[ch[ j ]]])>0)
                {
                    .infoPrint2(details, 2, "Marg onto sep %i: {%s}\n", ch[ j ], .colstr(seps[[ch[ j ]]]))
                    sp.pot            <- tableMargin(cqpotList[[ i ]], seps[[ch[ j ]]])
                    ##cat(sprintf("......is.na: sp.pot=%i\n", any(is.na(sp.pot))))
                    cqpotList[[ch[ j ]]]  <- tableOp2(cqpotList[[ch[ j ]]], sp.pot, `*`)
                    .infoPrint(details, 4, { cat("Marginal:\n"); print (sp.pot) })
                }
            }
        }
    }
    
    attr(cqpotList, "pEvidence") <- normConst
    cqpotList
}




