## ##############################################################
##
#' @title Update components of Bayesian network
#'
#' @description Update components of Bayesian network.
#'
#' @name cpt-update
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' 
##
## ##############################################################
#'
#' @param object A `grain` object.
#' @param value A named list, see examples below.
#' 
#' @seealso \code{\link{grain}}, \code{\link[gRbase]{propagate}},
#'     \code{\link[gRbase]{triangulate}}, \code{\link[gRbase]{rip}},
#'     \code{\link[gRbase]{junctionTree}}
#' @references Søren Højsgaard (2012). Graphical Independence
#'     Networks with the gRain Package for R. Journal of Statistical
#'     Software, 46(10), 1-26.
#'     \url{http://www.jstatsoft.org/v46/i10/}.
#' 
#' @keywords utilities models

#' @examples
#' ## See the wet grass example at
#' ## https://en.wikipedia.org/wiki/Bayesian_network
#' 
#' yn <- c("yes", "no")
#' p.R <- cptable(~R, values=c(.2, .8), levels=yn)
#' p.S_R <- cptable(~S:R, values=c(.01, .99, .4, .6), levels=yn)
#' p.G_SR <- cptable(~G:S:R, values=c(.99, .01, .8, .2, .9, .1, 0, 1), levels=yn)
#' 
#' x <- compileCPT(p.R, p.S_R, p.G_SR)
#' x
#' wet.bn <- grain(x)
#' 
#' getgrain(wet.bn, "cpt")
#' getgrain(wet.bn, "cpt")$R
#' getgrain(wet.bn, "cpt")$S
#'
#' # Now update some cpt's
#' wet.bn2 <- setCPT(wet.bn, list(R=c(.3, .7), S=c(.1, .9, .7, .3)))
#' 
#' getgrain(wet.bn2, "cpt")$R
#' getgrain(wet.bn2, "cpt")$S
#' 

#' @export 
#' @rdname cpt-update
setCPT <- function(object, value){
    UseMethod("setCPT")
}

#' @export 
#' @rdname cpt-update
setCPT.cpt_grain <- function(object, value){
    if (!.is.named.list(value))
        stop("'value' must be a named list")

    vn <- names(getgrain(object, "cpt"))
    nn <- names(value)
    if (any((id <- is.na(match(nn, vn)))))
        stop("variable(s) ", toString(nn[id]), " not in network")   
    for (i in seq_along(nn)){
        v <- nn[i]
        z <- value[[i]]
        if (length(z) != length(getgrain(object, "cpt")[[v]]))
            stop("replacement value not correct length")                    
        object$cptlist[[v]][] <- z               
    }
    isCompiled(object) <- FALSE
    isPropagated(object) <- FALSE
    object
}

.is.named.list <- function(x){
    inherits(x, "list") && !is.null(names(x))
}


## cpt-update
## "setcpt<-" <- function(object, value){
##     UseMethod("setcpt<-")
## }

## #' @rdname cpt-update
## "setcpt<-.grain" <- function(object, value){
##     setCPT(object, value)
## }



