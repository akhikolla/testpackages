#######################################################################
#'
#' @title Return the dimension of a log-linear model
#' 
#' @description Return the dimension of a log-linear model given by
#'     the generating class 'glist'. If the model is decomposable and
#'     adjusted dimension can be found.
#'
#' @name loglin-dim
#'
#######################################################################
#'
#' @aliases dim_loglin dim_loglin_decomp
#' 
#' @param glist Generating class (a list) for a log-linear model. See
#'     'details' below.
#' @param tableinfo Specification of the levels of the variables. See
#'     'details' below.
#' @param adjust Should model dimension be adjusted for sparsity of
#'     data (only available for decomposable models)
#'
#' @details
#' 
#' \code{glist} can be either a list of vectors with variable names or a list
#' of vectors of variable indices.
#' 
#' \code{tableinfo} can be one of three different things.
#' 
#' 1) A contingency table (a \code{table}).
#' 
#' 2) A list with the names of the variables and their levels (such as one
#' would get if calling \code{dimnames} on a \code{table}).
#' 
#' 3) A vector with the levels. If \code{glist} is a list of vectors with
#' variable names, then the entries of the vector \code{tableinfo} must be
#' named.
#' 
#' If the model is decomposable it \code{dim_loglin_decomp} is to be preferred over
#' \code{dim_loglin} as the former is much faster.
#' 
#' Setting \code{adjust=TRUE} will force \code{dim_loglin_decomp} to calculated a
#' dimension which is adjusted for sparsity of data. For this to work,
#' \code{tableinfo} *MUST* be a table.
#' 
#' @return A numeric.
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{dmod}}, \code{\link{glm}}, \code{\link[MASS]{loglm}}
#' @keywords models
#'
#' @examples
#' ## glist contains variable names and tableinfo is a named vector:
#' dim_loglin(list(c("a", "b"), c("b", "c")), c(a=4, b=7, c=6))
#' 
#' ## glist contains variable names and tableinfo is not named:
#' dim_loglin(list(c(1, 2), c(2, 3)), c(4, 7, 6))
#' 
#' ## For decomposable models:
#' dim_loglin_decomp(list(c("a", "b"), c("b", "c")), c(a=4, b=7, c=6),adjust=FALSE)
#' 

#' @export
#' @rdname loglin-dim
dim_loglin <- function(glist, tableinfo){

  if (is.table(tableinfo)){
    dtab <- c(lapply(dimnames(table), length), recursive=TRUE)
  } else {
    if (is.list(tableinfo)){
      dtab <- c(lapply(tableinfo, length), recursive=TRUE)
    } else {
      dtab <- tableinfo
    }
  }
  .dim_loglin(glist, dtab)
}

## Find dimension of decomposable model
## (with or without dimension adjustment for sparsity)
##
## 'cliq', 'seps' are cliques and separators (can be found from rip() function)
## 'table' can be either an array or a vector of levels with named components
##

#' @export
#' @rdname loglin-dim
dim_loglin_decomp <- function(glist, tableinfo, adjust=TRUE){
    rr <- ripMAT(.glist2adjMAT(glist))  ## FIXME glist2adjMAT should go...
    .dim_loglin_decomp(rr$cliques, rr$separators, tableinfo=tableinfo, adjust=adjust)
}

### dot-functions below here

.dim_loglin <- function(glist, dtab){

    .subsets <- function(x) {
        y <- list(vector(mode(x), length = 0))
        for (ii in seq_along(x)) {
            y <- c(y, lapply(y, "c", x[ii]))
        }
        y[-1]
    }
    
    max.g.size <- max(unlistPrim(lapply(glist, length)))
    ##cat("max.g.size:", max.g.size, "\n")

    if (max.g.size < 10)
    {
        ##print(lapply(glist, .subsets))
        zz    <- unlist(lapply(glist, .subsets), recursive=FALSE)
        unzz  <- unique.default(zz)
        ##print(unzz)
        
    }
    else
    {
        unzz <- .subsets(glist[[1]])
        base.idx  <- unlistPrim(lapply(unzz, function(terms) sum(2^(terms - 1)) ))
        if (length(glist)>1){
            for (ii in 2:length(glist)){
                tmp      <- .subsets(glist[[ii]])
                tmp.idx  <- unlistPrim(lapply(tmp, function(terms) sum(2^(terms - 1)) ))
                unzz     <- c(unzz,tmp[!(tmp.idx %in% base.idx)])
                base.idx <- c(base.idx, tmp.idx[!(tmp.idx %in% base.idx)])
            }
        }
    }
        
    ans <- 0
    for (jj in 1:length(unzz)){
        inc <- prod(dtab[unzz[[jj]]]-1)
        ans <- ans + inc
    }
    ans
}



.dim_loglin_decomp <- function(cliq, seps, tableinfo, adjust=TRUE){
    
    if (!adjust){
        if (is.array(tableinfo))
            vlev <- c(lapply(dimnames(tableinfo), length), recursive=TRUE)
        else
            vlev <- tableinfo
        
        ## Without adjustment of dimension for sparsity
        npar <- prod(vlev[cliq[[1]]])-1
        if (length(cliq) > 1){
            for (ii in 2:length(cliq)){
                dimC <- prod(vlev[cliq[[ii]]])-1
                dimS <- prod(vlev[seps[[ii]]])-1
                npar <- npar +  dimC - dimS
            }
        }
    } else {
    ## cat("With adjustment of dimension for sparsity\n")
        if (!is.array(tableinfo))
            stop("Model dimension adjusted for sparsity requires an array\n")
        tm1  <- tabMarg(tableinfo, cliq[[1]])
        npar <- sum(1 * (tm1 > 0)) - 1
        if (length(cliq) > 1){
            for (ii in 2:length(cliq)){
                tm1  <- tabMarg(tableinfo, cliq[[ii]])
                tm2  <- tabMarg(tm1, seps[[ii]])
                dimC <- sum(1 * (tm1 > 0)) - 1
                dimS <- sum(1 * (tm2 > 0)) - 1
                npar <- npar + dimC - dimS
            }
        }
    }
    return(npar)
    
    
    ## Extract from loglin().
    ## Just used for checking purposes
    ##
    .getDimkh <- function(glist, dtab){
        
        subsets <- function(x) {
            y <- list(vector(mode(x), length = 0))
            for (i in seq_along(x)) {
                y <- c(y, lapply(y, "c", x[i]))
            }
            y[-1]
        }
        
        nvar <- length(dtab)
        df <- rep.int(0, 2^nvar)
        for (k in seq_along(glist)) {
            terms <- subsets(glist[[k]])
            for (j in seq_along(terms)){
                df[sum(2^(terms[[j]] - 1))] <- prod(dtab[terms[[j]]] -  1)
            }
        }
        sum(df)
    }
}

