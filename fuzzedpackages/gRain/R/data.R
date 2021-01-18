#' Chest clinic example
#'
#' Conditional probability tables for the chest clinic example.
#' 
#' @name chest
#' @docType data
#'
#' @keywords datasets
#' @usage data(chest_cpt)
#'
#' @examples
#'
#' ## 'data' generated with the following code fragment
#' yn   <- c("yes", "no")
#' a    <- cptable(~asia, values=c(1,99),levels=yn)
#' t.a  <- cptable(~tub|asia, values=c(5,95,1,99),levels=yn)
#' s    <- cptable(~smoke, values=c(5,5), levels=yn)
#' l.s  <- cptable(~lung|smoke, values=c(1,9,1,99), levels=yn)
#' b.s  <- cptable(~bronc|smoke, values=c(6,4,3,7), levels=yn)
#' e.lt <- cptable(~either|lung:tub,values=c(1,0,1,0,1,0,0,1),levels=yn)
#' x.e  <- cptable(~xray|either, values=c(98,2,5,95), levels=yn)
#' d.be <- cptable(~dysp|bronc:either, values=c(9,1,7,3,8,2,1,9), levels=yn)
#'
#' grain(compileCPT(a, t.a, s, l.s, b.s, e.lt, x.e, d.be))
#' 
#' # 'data' generated from
#' # chest_cpt <- list(a, t.a, s, l.s, b.s, e.lt, x.e, d.be)
#' 
#' data(chest_cpt)
#' 
"chest_cpt"



#' Wet grass example
#'
#' Conditional probability tables for the wet grass example.
#' 
#' @name grass
#' @docType data
#'
#' @keywords datasets
#' @usage data(grass_cpt)
#'
#' @examples
#' 
#' ## 'data' generated with the following code fragment
#' yn <- c("yes", "no")
#' p.R    <- cptable(~R, values=c(.2, .8), levels=yn)
#' p.S_R  <- cptable(~S:R, values=c(.01, .99, .4, .6), levels=yn)
#' p.G_SR <- cptable(~G:S:R, values=c(.99, .01, .8, .2, .9, .1, 0, 1), levels=yn)
#' 
#' grain(compileCPT(p.R, p.S_R, p.G_SR))
#' 
#' # 'data' generated from
#' # grass_cpt <- list(p.R, p.S_R, p.G_SR)
#' 
#' data(grass_cpt)
#' 
"grass_cpt"
