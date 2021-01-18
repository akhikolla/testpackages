
#####################################################################
#' @title General functions related to iModels
#' @description General functions related to iModels
#' @name imodel-general
#####################################################################
#'
#' @param object,fit,x An \code{iModel} object.
#' @param scale Unused (and irrelevant for these models)
#' @param k Weight of the degrees of freedom in the AIC formula
#' @param ... Currently unused.

#' @rdname imodel-general
logLik.iModel <- function(object, ...){
    val <- fitinfo(object)$logL
    attr(val, "df") <- unname(fitinfo(object)$dimension["mod.dim"] )
    attr(val, "nobs") <- sum(datainfo(object)$data)
    class(val) <- "logLik"
    val
}

#' @rdname imodel-general
extractAIC.iModel <- function(fit, scale, k = 2, ...){
    unname(c(fitinfo(fit)$dimension["mod.dim"],
             -2 * fitinfo(fit)$logL + k * fitinfo(fit)$dimension["mod.dim"]))
}

#' @rdname imodel-general
summary.iModel <- function(object, ...){
  glist <- getmi(object, "glist")
  isg   <- getmi(object, "isGraphical")
  isd   <- getmi(object, "isDecomposable")

  cq    <- getCliques(ugList(glist))# $maxCliques
  ans   <- structure(list(glist=glist, isGraphical=isg, isDecomposable=isd,
                          cliques=cq),
                     class="iModelsummary")
  ans
}

#' @rdname imodel-general
print.iModelsummary <- function(x,...){
    cat(sprintf("is graphical=%s; is decomposable=%s\n",
                x$isGraphical, x$isDecomposable))
    cat("generators (glist):\n")
    str(.glist(x), give.head=FALSE, comp.str=" ", no.list=TRUE)
    ##cat("EXPERIMENTAL: components: ", names(x),"\n")
    invisible(x)
}

.extractFIT <- function(object,...){
  c(object[[1]], object$df)
}

.glist2formula <- function (f) {
    if (inherits(f, "formula"))
        return(f)
    ans <- try(as.formula(paste("~", paste(unlist(lapply(f, paste, collapse = "*")),
                                           collapse = "+")), .GlobalEnv),silent=TRUE)
    if (inherits(ans, "try-error"))
        stop("Unable to create formula from list. \nCould be due to white space, strange characters etc. in variable names\n")
    ans
}

#' importFrom stats formula terms

#' @export
#' @rdname imodel-general
formula.iModel <- function(x,...){
    .glist2formula(terms(x))
}

#' @export
#' @rdname imodel-general
terms.iModel <- function(x, ...){
    modelinfo(x)$glist
}


#' @rdname imodel-general
isGraphical.dModel <- function(x){
    isGraphical(terms(x))
}

#' @rdname imodel-general
isDecomposable.dModel <- function(x){
    isDecomposable(terms(x))
}

#' @rdname imodel-general
modelProperties <- function(object){
    UseMethod("modelProperties")
}

#' @rdname imodel-general                                 
modelProperties.dModel <- function(object){
    x <- terms( object )
    vn <- unique(unlist(x))
    amat <- .glist2adjMAT(x, vn = vn)  ## FIXME glist2adjMAT
    cliq <- maxCliqueMAT(amat)[[1]]
    ##isg <- all(unlist(lapply(cliq, function(cq) isin(x, cq))))  ## FIXME isin
    isg <- all(unlist(lapply(cliq, function(cq) is_inset(cq, x))))  ## FIXME isin
    isd <- if (isg) {
               length(mcsMAT(amat)) > 0
           }
           else FALSE
    
    c(isGraphical=isg, isDecomposable=isd)
}


datainfo <- function(object){
    object$datainfo
}

fitinfo <- function(object){
    object$fitinfo
}

modelinfo <- function(object){
    object$modelinfo
}
