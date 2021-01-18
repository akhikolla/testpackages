##FIXME: Delete df from loglin() output.

##########################################################
##
#' @title Discrete interaction model (log-linear model)
#'
#' @description Specification of log--linear (graphical) model. The
#'     'd' in the name \code{dmod} refers to that it is a (graphical)
#'     model for 'd'iscrete variables
#' 
#' @name imodel-dmod
##
##########################################################

#' @details The independence model can be specified as \code{~.^1} and
#'     \code{~.^.} specifies the saturated model.  Setting
#'     e.g. \code{interactions=3} implies that there will be at most
#'     three factor interactions in the model.
#' 
#' Data can be specified as a table of counts or as a dataframe. If
#' data is a dataframe then it will be converted to a table (using
#' \code{xtabs()}). This means that if the dataframe contains numeric
#' values then the you can get a very sparse and high dimensional
#' table. When a dataframe contains numeric values it may be
#' worthwhile to discretize data using the \code{cut()} function.
#' 
#' The \code{marginal} argument can be used for specifying the
#' independence or saturated models for only a subset of the
#' variables. When \code{marginal} is given the corresponding marginal
#' table of data is formed and used in the analysis (notice that this
#' is different from the behaviour of \code{loglin()} which uses the
#' full table.
#' 
#' The \code{triangulate()} method for discrete models (dModel
#' objects) will for a model look at the dependence graph for the
#' model.
#' 
#' @aliases dmod print.dModel fitted.dModel residuals.dModel 
#'     triangulate.dModel
#' @param formula Model specification in one of the following forms: 1) a
#'     right-hand sided formula, 2) as a list of generators, 3) an undirected
#'     graph (represented either as a graphNEL object or as an adjacency
#'     matrix).  Notice that there are certain model specification shortcuts,
#'     see Section 'details' below.
#' @param data Either a table or a dataframe. In the latter case, the dataframe
#'     will be coerced to a table. See 'details' below.
#' @param interactions A number given the highest order interactions in the
#'     model, see Section 'details' below.
#' @param marginal Should only a subset of the variables be used in connection
#'     with the model specification shortcuts
#' @param fit Should the model be fitted.
#' @param details Control the amount of output; for debugging purposes.
#' @param ... Additional arguments; currently no used.
#'
#' @return An object of class \code{dModel}.
#'
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#'
#' @seealso \code{\link{cmod}}, \code{\link{mmod}}
#' @keywords models
#' @examples
#'  
#' ## Graphical log-linear model
#' data(reinis)
#' dm1 <- dmod(~ .^., reinis)
#' dm2 <- backward(dm1, k=2)
#' dm3 <- backward(dm1, k=2, fixin=list(c("family", "phys", "systol")))
#' ## At most 3-factor interactions
#' dm1<-dmod(~ .^., data=reinis, interactions=3)


#' @export dmod
dmod <- function(formula, data, marginal=NULL, interactions=NULL, fit=TRUE, details=0, ...){

    if (!inherits(data, c("data.frame", "table")))
      stop("data must be either a dataframe or a table \n")
  
    if (inherits(data, "data.frame")){
        is_data.frame <- TRUE
        varNames <- names(data)
    } else {
        is_data.frame <- FALSE
        varNames <- names(dimnames(data))
    }
    
    if (length(marginal) > 0){
        zzz  <- unlist(lapply(marginal, pmatch, varNames))
        zzz  <- zzz[!is.na(zzz)]
        marginal <- varNames[zzz]
    }
    
    ans <- parse_gm_formula(formula, varNames, marginal, interactions)
    
    
    if (is_data.frame){
        data <- xtabs(~., data=data[ans$varNames])
    } else {
        if (length(ans$varNames) != length(varNames)){
            ## FIXME: Looks strange: as.table(tabMarg(data, ans$varNames))
            data <- as.table(tabMarg(data, ans$varNames))
        }
    }
    varNames     <- names(dimnames(data))
    
    res <- list(
        modelinfo = NULL,
        varNames  = varNames,
        datainfo  = list(data=data),
        fitinfo   = NULL,
        isFitted  = FALSE
    )
    
    upd <- .dModel_finalize(ans$glist, varNames) ## NOTE not .glist
    res$modelinfo  <- upd    

    class(res) <- c("dModel", "iModel")
    if (fit) fit(res) else res
}


.dModel_finalize<- function(glist, varNames){
    list(glist      = glist,
         glistNUM   = .glistNUM(glist, varNames),
         properties = isGSD_glist(glist))
}





#' @export
fitted.dModel <- function(object,...){
    if ( object$isFitted ){
        object$fitinfo$fit
    } else {
        cl <- object$call
        cl$fit <- TRUE
        eval( cl )$fitinfo$fit
    }
}


#' @export
fit.dModel <- function(object, engine="loglin", print=FALSE, ...){

  ## FIXME: At some point we should allow for data in the form of a dataframe
  ##
  switch(engine,
         "loglin"={llfit <- loglin(getmi(object, "data"), getmi(object, "glist"),
                                   fit=TRUE,  print=print, eps=0.10, ...)
                                   names(llfit)[1] <- 'dev'  ## loglin() calls this slot 'lrt'
               })

    vn <- getmi(object, "varNames")
    
    ## Calculate df's and adjusted df's. Requires data, glist
    glistNUM <- .glistNUM(getmi(object, "glist"), vn)
    
    sat.dim.unadj   <- prod(dim(getmi(object, "data"))) - 1
    sat.dim.adj     <- sum(getmi(object, "data") > 0) - 1
    ind.dim.unadj   <- sum(dim(getmi(object, "data")) - 1)
    
    if (getmi(object, "isDecomposable")){
        rr <- ripMAT(.glist2adjMAT(getmi(object, "glist")))
        dim.adj   <- .dim_loglin_decomp(rr$cliques, rr$separators,
                                   tableinfo=getmi(object, "data"), adjust=TRUE)
        dim.unadj <- .dim_loglin_decomp(rr$cliques, rr$separators,
                                   tableinfo=getmi(object, "data"), adjust=FALSE)
    } else {
        dim.adj   <- NA
        dim.unadj <- .dim_loglin(glistNUM, dim(getmi(object, "data")))
    }
    
    df.adj       <- sat.dim.adj   - dim.adj
    df.unadj     <- sat.dim.unadj - dim.unadj

    ## FIXME: SILLY to call loglin to fit independence model, but it is faster than
    ## a simple R implementation.
    #indep.model  <- loglin(getmi(object, "data"), as.list(vn), iter=1, print=FALSE)
    #ideviance    <-  -(llfit$dev - indep.model$lrt)
    #idf          <-  -(llfit$df - indep.model$df)

    indep.stat  <- loglin(getmi(object, "data"), as.list(vn), iter=1, print=FALSE)[c("lrt", "df")]
    ideviance    <-  -(llfit$dev - indep.stat$lrt)
    idf          <-  -(llfit$df - indep.stat$df)
    
    ii <- getmi(object, "data") * llfit$fit > 0
    llfit$logL      <- sum(getmi(object, "data")[ii] * log(llfit$fit[ii] / sum(llfit$fit)))
    llfit$ideviance <- ideviance

    extra1 <- list(dim.unadj = dim.unadj,
                   dim.adj   = dim.adj,
                   df.unadj  = df.unadj,
                   df.adj    = df.adj,
                   idf       = idf,
                   ideviance = ideviance)
    
    dimension <- c(mod.dim=dim.unadj, sat.dim=sat.dim.unadj,
                   i.dim=ind.dim.unadj, df=df.unadj, idf=idf,
                   mod.dim.adj = dim.adj,
                   sat.dim.adj = sat.dim.adj,
                   df.adj      = df.adj  )
    
    llfit$df        <- NULL
    llfit$fit <- NULL
    llfit$aic       <- -2 * llfit$logL + 2 * dimension['mod.dim']
    llfit$bic       <- -2 * llfit$logL + log(sum(getmi(object, "data"))) * dimension['mod.dim']
    
    ## Calculate warning codes for sparsity
    sparseinfo <- .dModel_sparsity (llfit,
                                    getmi(object, "isDecomposable"),
                                    getmi(object, "glist"),
                                    getmi(object, "data"))
    
    fitinfo <- llfit
    fitinfo$sparseinfo <- sparseinfo
    fitinfo$dimension  <- dimension
    
    object$fitinfo  <- fitinfo
    object$isFitted <- TRUE
    class(object)   <- c("dModel", "iModel")
    object
}


.dModel_sparsity <- function(llfit, isDecomposable, glist, data){

    ## Calculate warning codes for sparsity
  ##
  fc       <- as.numeric(llfit$fit)
  all.gt.0 <- sum( fc < 0.00001 ) == 0
  all.gt.5 <- sum( fc < 5 ) == 0

  frac.gt.0 <- sum(fc>0)/length(fc)
  frac.gt.5 <- sum(fc>5)/length(fc)

  if (!isDecomposable)
    {
      ## cat ("Model is not decompsable...\n")
      if (all.gt.0){
        df.ok        <- TRUE
        sparse.df.ok <- FALSE
        all.gt.5 <- sum( fc < 5 ) == 0
        if (!all.gt.5){
          #cat("Warning: table is sparse and asymptotic chi2 distribution is questionable (2)\n")
          chi2.ok <- FALSE
        } else {
          chi2.ok <- TRUE
        }
      } else {
        #cat("degrees of freedom are unadjusted and can not be trusted (1)\n")
        #cat("Warning: asymptotic chi2 distribution is questionable (2)\n")
        df.ok        <- FALSE
        sparse.df.ok <- FALSE
        chi2.ok      <- FALSE
      }
    }
  else
    {
      ## cat("Model is decompsable...\n")
      if (all.gt.0){
        df.ok         <- TRUE
        sparse.df.ok  <- TRUE
      } else {
        df.ok         <- FALSE
        sparse.df.ok  <- TRUE
      }

      all.gt.5 <- 1 - sum( fc < 5 ) > 0
      if (all.gt.5){
        ## It is a dense table
        chi2.ok <- TRUE
      } else {
        ## It is not a dense table, but the table may be dense in the cliques
        ## Then find, in each clique, those marginal cells with positive counts. This leads to the df adjustment.
        ## For each marginal cell with positive counts, check if counts are > 5. If so, the chi2 is ok.
        sparsecode3 <- rep.int(0, length(glist))
        for (ii in 1:length(glist)){
          tmC <- tabMarg(data, glist[[ii]])
          tm0 <- tmC > 0
          tm5 <- tmC[tm0] > 5
          if (length(tm5) == length(tm0)){
            sparsecode3[ii] <- 1
          }
        }
        all.gt.5.when.gt.0 <- sum(sparsecode3)==length(glist)

        if (all.gt.5.when.gt.0){
          #cat("Warning: table is sparse and degrees of freedom have been adjusted to reflect sparsity of table (3)\n")
          chi2.ok <- TRUE
        } else {
          #cat("Warning: table is sparse and asymptotic chi2 distribution is questionable (4)\n")
          chi2.ok <- FALSE
        }
      }
    }

  sparseinfo <-
    c(chi2.ok=chi2.ok, df.ok=df.ok, sparse.df.ok=sparse.df.ok)

  sparseinfo
}


#' @method residuals "dModel"
#' @export
residuals.dModel <-
    function (object, type = c("deviance", "pearson", "response"),
              ...)
{
    type <- match.arg(type)
    y  <- getmi(object, "data")
    mu <- fitted(object)
    res <- y - mu
    nz <- mu > 0
    y <- y[nz]
    mu <- mu[nz]
    res[nz] <- switch(type, deviance = sign(y - mu) * sqrt(2 *
        abs(y * log((y + (y == 0))/mu) - (y - mu))), pearson = (y -
        mu)/sqrt(mu), response = y - mu)
    res
}



