#' @title Create conditional probability tables (CPTs)
#' 
#' @description Creates conditional probability tables of the form
#'     p(v|pa(v)).
#'
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' 
#' @param vpar Specifications of the names in P(v|pa1,...pak). See
#'     section 'details' for information about the form of the
#'     argument.
#' @param values Probabilities; recycled if necessary. Regarding the
#'     order, please see section 'details' and the examples.
#' @param normalize See 'details' below.
#' @param smooth See 'details' below.
#' @param levels See 'details' below.
#'
#' @details
#' 
#' If \code{normalize=TRUE} then the probabilities are normalized to sum to one
#' for each configuration of the parents.
#' 
#' If \code{smooth} is non--zero then zero entries of \code{values} are
#' replaced with \code{smooth} __before__ normalization takes place.
#' 
#' Regarding the form of the argument \code{vpar}: To specify \eqn{P(a|b,c)}
#' one may write \code{~a|b:c}, \code{~a:b:c}, \code{~a|b+c}, \code{~a+b+c} or
#' \code{c("a","b","c")}. Internally, the last form is used. Notice that the
#' \code{+} and \code{:} operator is used as a separator only. The order of the
#' variables IS important so the operators DO NOT commute.
#' 
#' If \code{a} has levels \code{a1,a2} and likewise for \code{b} and \code{c}
#' then the order of \code{values} corresponds to the configurations
#' \code{(a1,b1,c1)}, \code{(a2,b1,c1)} \code{(a1,b2,c1)}, \code{(a2,b2,c1)}
#' etc. That is, the first variable varies fastest.  Hence the first two
#' elements in \code{values} will be the conditional probabilities of \code{a}
#' given \code{b=b1, c=c1}.
#' 
#' @return A \code{cptable} object (a numeric vector with various attributes).
#' 
#' @seealso \code{\link{andtable}}, \code{\link{ortable}},
#'     \code{\link{extractCPT}}, \code{\link{compileCPT}},
#'     \code{\link{extractPOT}}, \code{\link{compilePOT}},
#'     \code{\link{grain}}, \code{\link[gRbase]{parray}}
#' @references Søren Højsgaard (2012). Graphical Independence
#'     Networks with the gRain Package for R. Journal of Statistical
#'     Software, 46(10), 1-26.
#'     \url{http://www.jstatsoft.org/v46/i10/}.
#' @keywords models
#' @examples
#' 
#' 
#' ## See the wet grass example at
#' ## https://en.wikipedia.org/wiki/Bayesian_network
#' 
#' yn <- c("yes", "no")
#' p.R    <- cptable(~R, values=c(.2, .8), levels=yn)
#' p.S_R  <- cptable(~S:R, values=c(.01, .99, .4, .6), levels=yn)
#' p.G_SR <- cptable(~G:S:R, values=c(.99, .01, .8, .2, .9, .1, 0, 1), levels=yn)
#'
#' # or
#' ssp <- list(R=yn, S=yn, G=yn) # state space
#' p.R    <- cptable(~R, values=c(.2, .8), levels=ssp)
#' p.S_R  <- cptable(~S:R, values=c(.01, .99, .4, .6), levels=ssp)
#' p.G_SR <- cptable(~G:S:R, values=c(.99, .01, .8, .2, .9, .1, 0, 1), levels=ssp)
#'
#' # components above are "intermediate representations" and are turned into arrays with
#' wet.cpt <- compileCPT(p.R, p.S_R, p.G_SR)
#' wet.cpt
#' wet.cpt$S # etc
#'
#' # A Bayesian network is created with:
#' wet.bn <- grain(wet.cpt)
#' 
#' # Can also create arrays directly
#' \dontrun{
#' ssp <- list(R=yn, S=yn, G=yn) # state space
#' p.R    <- c(.2, .8)
#' p.S_R  <- c(.01, .99, .4, .6)
#' p.G_SR <- c(.99, .01, .8, .2, .9, .1, 0, 1)
#' dim(p.R) <- 2
#' dimnames(p.R) <- ssp["R"]
#' dim(p.S_R) <- c(2, 2)
#' dimnames(p.S_R) <- ssp[c("S", "R")]
#' dim(p.G_SR) <- c(2, 2, 2)
#' dimnames(p.G_SR) <- ssp[c("G", "S", "R")]
#'
#' # Arrays can be created (easier?) with parray() from gRbase
#' p.R    <- parray("R", levels=ssp, values=c(.2, .8))
#' p.S_R  <- parray(c("S", "R"), levels = ssp, values=c(.01, .99, .4, .6))
#' p.G_SR <- parray(~ G:S:R, levels = ssp, values=c(.99, .01, .8, .2, .9, .1, 0, 1))
#' }



#' @export
cptable <- function(vpar, levels=NULL, values=NULL, normalize=TRUE,  smooth=0 ){
    vpa  <- c(.formula2char(vpar))        
    if (is.list(levels)){
        v <- vpa[1]
        if (!(v %in% names(levels)))
            stop(paste0("Name ", v, " is not in the 'levels' list\n"))
        levels <- levels[[v]]
    }
    ##str(list(vpa=vpa, xlevels=levels))
    ## if (is.null(values))
    ##     values <- rep(1.0, length(levels))
    out  <- values
    attributes(out) <-
        list(vpa=vpa, normalize=normalize,
             smooth=smooth, levels=levels)
    class(out) <- "cptable_class"
    out
}


## #' @rdname cptable
## cptab <- cptable


## NORMAL

cdist <- function(vpar, parm=list()){
    vpa  <- c(.formula2char(vpar))        
    out <- parm
    names(out) <- c("intercept", "slope", "sd")
    attr(out, "vpa") <- vpa
    class(out) <- "cdist"
    out
}


#' @export
print.cptable_class <- function(x, ...){
    ## "print.cptable\n" %>% cat
    v <- c(x)
    dim(v) <- c(length(attr(x,"levels")), length(v) / length(attr(x, "levels")))
    rownames(v) <- attr(x, "levels")
    cat(sprintf("{v, pa(v)} :\n"))
    str(attr(x, "vpa"))
    print(v)    
    ##str(attributes(x))
    invisible(x)
}


summary.cptable_class <- function(object, ...){
    print(object)
    str(attributes(object))
    invisible(object)
}

varNames.cptable_class <- function(x){
    ##x$vpa
    attr(x, "vpa")
}

valueLabels.cptable_class <- function(x){
    out <- list(attr(x, "levels"))
    nam <- attr(x, "vpa")
    names(out) <- attr(x, "vpa")[1] #x$vpa[1]
    out
}



