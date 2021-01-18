#' @title Query a network
#' 
#' @description Query an independence network, i.e. obtain the
#'     conditional distribution of a set of variables - possibly (and
#'     typically) given finding (evidence) on other variables.
#' 
#' @name querygrain
#' 
#' @aliases querygrain querygrain.grain qgrain
#'
#' @param object A `grain` object.
#' @param nodes A vector of nodes; those nodes for which the
#'     (conditional) distribution is requested.
#' @param evidence An alternative way of specifying findings
#'     (evidence), see examples below.
#' @param exclude If \code{TRUE} then nodes on which evidence is given
#'     will be excluded from \code{nodes} (see above).
#' @param normalize Should the results be normalized to sum to one.
#' @param type Valid choices are \code{"marginal"} which gives the
#'     marginal distribution for each node in \code{nodes};
#'     \code{"joint"} which gives the joint distribution for
#'     \code{nodes} and \code{"conditional"} which gives the
#'     conditional distribution for the first variable in \code{nodes}
#'     given the other variables in \code{nodes}.
#' @param result If "data.frame" the result is returned as a data
#'     frame (or possibly as a list of dataframes).
#' @param details Debugging information
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
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{setEvidence}}, \code{\link{getEvidence}},
#'     \code{\link{retractEvidence}}, \code{\link{pEvidence}}
#' @references Søren Højsgaard (2012). Graphical Independence Networks
#'     with the gRain Package for R. Journal of Statistical Software,
#'     46(10), 1-26.  \url{http://www.jstatsoft.org/v46/i10/}.
#' @keywords models utilities
#' @examples
#' 
#' testfile <- system.file("huginex", "chest_clinic.net", package = "gRain")
#' chest <- loadHuginNet(testfile, details=0)
#' qb <- querygrain(chest)
#' qb
#' 
#' lapply(qb, as.numeric) # Safe
#' sapply(qb, as.numeric) # Risky
#' 
#' @export
querygrain <- function(object, nodes=nodeNames(object), type="marginal",
                       evidence=NULL, exclude=TRUE, normalize=TRUE,
                       result="array", details=0)
{
  UseMethod("querygrain")
}

#' @export 
qgrain <- querygrain

#' @export
querygrain.grain <- function(object, nodes = nodeNames(object), type = "marginal",
                             evidence=NULL, exclude=TRUE, normalize=TRUE,
                             result="array", details=0){

    if (!is.null(evidence)){
        if (details>=1) cat(" Inserting (additional) evidence\n")
        object <- setEvidence(object, evidence=evidence)
    }
    
    type <- match.arg(type, c("marginal", "joint", "conditional"))
    result <- match.arg(result, c("array","data.frame"))
    t0 <- proc.time()

    if (is.null(nodes)) return(invisible(NULL))
    else if (inherits(nodes, "formula")) nodes <- unlist(rhsf2list(nodes))

    if (!isCompiled(object)){ 
        if (details >= 1) cat("  Compiling (and propagating) model ...\n")
        object <- compile(object, propagate=TRUE)
    } else {
        if (!isPropagated(object)){
            if (details>=1) cat("  Propagating model...\n")
            object <- propagate(object)
        }
    }

    type = match.arg(type, choices=c("marginal", "joint", "conditional"))
    switch(type,
           "marginal"={
               out <- .nodeMarginal(object, nodes=nodes, exclude=exclude,
                                    details=details)
           },
           "joint"={
               out <- .nodeJoint(object, nodes=nodes, exclude=exclude,
                                 normalize=normalize, details=details)
           },
           "conditional"={
               qobject <- querygrain.grain(object, nodes=nodes,
                                           type="joint", exclude=exclude)
               if (length(names(dimnames(qobject))) > 1)
                   out <- tabDiv( qobject, tabMarg(qobject, nodes[-1]) )
               else
                   out <- qobject
           })

    if (result == "data.frame")
        out <- lapply(out, as.data.frame.table)

    if (object$control$timing)
        cat("Time: query", proc.time()-t0, "\n")
    out
}

.nodeJoint <- function(object, nodes=NULL, exclude=TRUE,
                       normalize=TRUE, details=1){

    nodes <- if (is.null(nodes)){rip(object)$nodes}
             else intersect(rip(object)$nodes, nodes)
    

    if (exclude)
        nodes <- setdiff(nodes, getEvidence(object)$nodes)

    cliq  <- rip(object)$cliques
    ##idxb <- sapply(cliq, function(cq) is_subsetof_(nodes, cq))
    idxb <- sapply(cliq, function(cq) is_subsetof(nodes, cq)) 

    if (any(idxb)){
        ##cat(".Calculating joint directly from clique\n")
        pt2   <- getgrain(object, "pot_equi")[[ which(idxb)[1] ]]
        value <- tableMargin(pt2, nodes)  ## FIXME tableMargin

        if (!normalize){
            zz <- value * pEvidence(object)
            value[] <- zz
        }

    } else {
        ##cat(".Calculating marginal brute force\n")
        nnodes <- length(nodes)
        dn     <- universe(object)$levels[nodes]
        value  <- tabNew(names(dn), dn) 

        nodes2 <- nodes[2:nnodes]
        dn2    <- universe(object)$levels[nodes2]
        gr2    <- as.matrix( expand.grid( dn2 ) )

        object  <- absorbEvidence( object )
        zz <- lapply(1:nrow(gr2), function(i){
            tmp <- setFinding(object, nodes=nodes2, states=gr2[i,])
            r   <- .nodeMarginal(tmp, nodes[1])
            v   <- r[[1]] * pEvidence(tmp)
            v
        })
        zz <- unlist(zz, use.names=FALSE)
        zz[is.na( zz )] <- 0
        if (normalize)
            zz <- zz / sum( zz )
        value[] <- zz
    }
    value
}

##idxb <- sapply(cliq, function(cq) subsetof(nodes, cq))
##tab   <- pot(object)$pot_equi[[ which(idxb)[1] ]]


.nodeMarginal <- function(object, nodes=NULL, exclude=TRUE, details=1){
    ##cat(".nodeMarginal\n")

    nodes <- if (is.null(nodes)){rip(object)$nodes}
             else intersect(rip(object)$nodes, nodes)
    
    if (exclude)
        nodes <- setdiff(nodes, getEvidence(object)$nodes)

    idx <- match(nodes, rip(object)$nodes)
    host.cq <- rip(object)$host[idx]

    if (length(nodes) > 0){
        out <- vector("list", length(nodes))
        names(out) <- nodes

        for (i in 1:length(nodes)){
            cvert  <- nodes[i]
            idx    <- host.cq[i]
            ## querygrain - .nodeMarginal: Calculations based on equipot
                                        #cpot   <- pot(object)$pot_equi[[ idx ]]
            cpot   <- getgrain(object, "pot_equi")[[ idx ]]
            ##mtab   <- tableMargin( cpot, cvert ) ## FIX tableMargin replaced
            mtab <- tabMarg(cpot, cvert)
            
            mtab   <- mtab / sum( mtab )
            out[[ i ]] <- mtab
        }
        return( out )
    }
}
