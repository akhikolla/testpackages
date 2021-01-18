###########################################################################
#'
#' @title Graphical Gaussian model
#' 
#' @description Specification of graphical Gaussian model. The 'c' in
#'     the name \code{cmod} refers to that it is a (graphical) model
#'     for 'c'ontinuous variables
#'
###########################################################################
#'
#' @details The independence model can be specified as \code{~.^1} and
#'     the saturated model as \code{~.^.}. The \code{marginal}
#'     argument can be used for specifying the independence or
#'     saturated models for only a subset of the variables.
#' @aliases extract_cmod_data
#' @param formula Model specification in one of the following forms:
#'     1) a right-hand sided formula, 2) as a list of generators.
#'     Notice that there are certain model specification shortcuts,
#'     see Section 'details' below.
#' 
#' @param data Data in one of the following forms: 1) A dataframe or
#'     2) a list with elements \code{cov} and \code{n.obs} (such as
#'     returned by the \code{cov.wt()} function.)
#'
#' @param marginal Should only a subset of the variables be used in
#'     connection with the model specification shortcuts.
#'
#' @param fit Should the model be fitted.
#'
#' @param details Control the amount of output; for debugging
#'     purposes.
#' @return An object of class \code{cModel} (a list)
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{dmod}}, \code{\link{mmod}},
#'     \code{\link{ggmfit}}
#' @keywords models
#' @examples
#' 
#' ## Graphical Gaussian model
#' data(carcass)
#' cm1 <- cmod(~ .^., data=carcass)
#' 
#' ## Stepwise selection based on BIC
#' cm2 <- backward(cm1, k=log(nrow(carcass)))
#' 
#' ## Stepwise selection with fixed edges
#' cm3 <- backward(cm1, k=log(nrow(carcass)),
#'  fixin=matrix(c("LeanMeat", "Meat11", "Meat12", "Meat13",
#'                 "LeanMeat", "Fat11", "Fat12", "Fat13"),
#'                  ncol=2))
#' 
#' @export cmod 

cmod <- function(formula, data, marginal=NULL, fit=TRUE, details=0){

    dd  <- extract_cmod_data(data)
    vn  <- colnames(dd$S)
    ans <- parse_gm_formula(formula, vn, marginal)
    ## Get varNames in the order matching to the data:
    vn <- vn[sort(match(ans$varNames, vn))]

    datainfo <- list(S=dd$S[vn, vn],
                     n.obs=dd$n.obs,
                     data=data)
    
    res <- list(modelinfo      = NULL,
                varNames       = vn,
                datainfo       = datainfo,
                fitinfo        = NULL,
                isFitted       = FALSE)
    
    upd   <- .cModel_finalize(ans$glist, vn)  ## NOTE not .glist here
    res$modelinfo <- upd  

    class(res) <- c("cModel", "iModel")
    if (fit) fit(res) else res
}


.cModel_finalize <- function(glist, varNames){

    amat  <- ugList(glist, result="matrix")
    glist <- maxCliqueMAT(amat)[[1]]
    isd   <- length(mcsMAT(amat)) > 0   
    
    list(glist       = glist,
         glistNUM    = .glistNUM(glist, varNames), 
         properties  = c(isg=TRUE, issd=isd))                
}




#' @export
fit.cModel <- function(object, engine="ggmfit", start=NULL, ...){

    fitfun <- if (identical(engine, "ggmfit")) ggmfit else ggmfitr

    ff <- fitfun(object$datainfo$S, n.obs=object$datainfo$n.obs,
                        glist=object$modelinfo$glist,
                        start=start, details=0,...)

    ## ideviance to independence model  
    idev  <-  ff$n.obs * (log(ff$detK) + sum(log(diag(ff$S))))  
    idim      <-  ff$nvar 
    sat.dim   <-  ((idim + 1) * idim) / 2
    dim.unadj <-  sat.dim - ff$df
    
    idf       <-  (dim.unadj - idim)
    logL.sat  <-  ff$logL + ff$dev/2
    
    aic       <-  -2 * ff$logL + 2 * dim.unadj
    bic       <-  -2 * ff$logL + log(ff$n.obs) * dim.unadj
    
    dimension <- c(mod.dim=dim.unadj, sat.dim=sat.dim, i.dim=idim, df=ff$df, idf=idf)
    
    ans   <- list(dev=ff$dev, ideviance=idev, logL.sat=logL.sat,
                  aic=aic, bic=bic,
                  dimension=dimension)
    
    ff$S <- ff$n.obs <- ff$dev <- ff$df <- NULL
    ans <- c(ff, ans)
    
    object$fitinfo  <- ans
    object$isFitted <- TRUE
    class(object)   <- c("cModel", "iModel")
    object
}


#' @export
extract_cmod_data <- function(data){
    if (inherits(data, "data.frame")){
        data <- cov.wt(data, method="ML")
    } else
        if (inherits(data, "list") && identical(names(data), c("cov", "center", "n.obs"))){
            ## OK
        } else
            stop("Can not proceed...")
    
            
    names(data)[1] <- "S"
    data
}











