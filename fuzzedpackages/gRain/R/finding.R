## #################################################################
##
## setFinding etc:
##
## Must live to ensure backward compatibility (it is used in the GMwR
## book, in the JSS paper and in many places on the web
##
## The functions call the corresponding evidence functions,
## setEvidence etc.
##
## #################################################################

#' @title Set, retrieve, and retract finding in Bayesian network.
#' 
#' @description Set, retrieve, and retract finding in Bayesian
#'     network.  NOTICE: The functions described here are kept only
#'     for backward compatibility; please use the corresponding
#'     evidence-functions in the future.
#' 
#' @name finding
#' 
#' @aliases setFinding retractFinding getFinding pFinding
#' @param object A "grain" object
#' @param nodes A vector of nodes
#' @param states A vector of states (of the nodes given by 'nodes')
#' @param flist An alternative way of specifying findings, see
#'     examples below.
#' @param propagate Should the network be propagated?
#' @note NOTICE: The functions described here are kept only for
#'     backward compatibility; please use the corresponding
#'     evidence-functions in the future:
#' 
#' \code{setEvidence()} is an improvement of \code{setFinding()} (and as such
#' \code{setFinding} is obsolete). Users are recommended to use
#' \code{setEvidence()} in the future.
#' 
#' \code{setEvidence()} allows to specification of "hard evidence" (specific
#' values for variables) and likelihood evidence (also known as virtual
#' evidence) for variables.
#' 
#' The syntax of \code{setEvidence()} may change in the future.
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' 
#' @seealso \code{\link{setEvidence}}, \code{\link{getEvidence}},
#' \code{\link{retractEvidence}}, \code{\link{pEvidence}},
#' \code{\link{querygrain}}
#'
#' @references Søren Højsgaard (2012). Graphical Independence Networks
#' with the gRain Package for R. Journal of Statistical Software, 46(10), 1-26.
#' \url{http://www.jstatsoft.org/v46/i10/}.
#'
#' @keywords models utilities
#' @examples
#' 
#' 
#' ## setFindings
#' yn <- c("yes", "no")
#' a    <- cptable(~asia, values=c(1,99),levels=yn)
#' t.a  <- cptable(~tub+asia, values=c(5,95,1,99),levels=yn)
#' s    <- cptable(~smoke, values=c(5,5), levels=yn)
#' l.s  <- cptable(~lung+smoke, values=c(1,9,1,99), levels=yn)
#' b.s  <- cptable(~bronc+smoke, values=c(6,4,3,7), levels=yn)
#' e.lt <- cptable(~either+lung+tub,values=c(1,0,1,0,1,0,0,1),levels=yn)
#' x.e  <- cptable(~xray+either, values=c(98,2,5,95), levels=yn)
#' d.be <- cptable(~dysp+bronc+either, values=c(9,1,7,3,8,2,1,9), levels=yn)
#' chest.cpt <- compileCPT(a, t.a, s, l.s, b.s, e.lt, x.e, d.be)
#' chest.bn <- grain(chest.cpt)
#' 
#' ## These two forms are equivalent
#' bn1 <- setFinding(chest.bn, nodes=c("chest", "xray"), states=c("yes", "yes"))
#' bn2 <- setFinding(chest.bn, flist=list(c("chest", "yes"), c("xray", "yes")))
#' 
#' getFinding(bn1)
#' getFinding(bn2)
#' 
#' pFinding(bn1)
#' pFinding(bn2)
#' 
#' bn1 <- retractFinding(bn1, nodes="asia")
#' bn2 <- retractFinding(bn2, nodes="asia")
#' 
#' getFinding(bn1)
#' getFinding(bn2)
#' 
#' pFinding(bn1)
#' pFinding(bn2)
#' 
#' 
#' @export setFinding
setFinding <- function(object, nodes=NULL, states=NULL, flist=NULL, propagate=TRUE){
    if (!is.null(flist)){
        flist2 <- do.call("rbind",flist)
        nodes   <- flist2[,1]
        states  <- flist2[,2]
    }
    setEvidence(object, nodes=nodes, states=states, propagate=propagate)
}

#' @export 
retractFinding <- retractEvidence

#' @export 
pFinding <- pEvidence

#' @export 
getFinding <- getEvidence

