## ====================================================================
## An object of class "covComp" has components kernels where some
## parameters are stored, and "top" parameters
##
## Doc generated in a classical (non-Roxygen) fashion.
## ====================================================================
setClass("covComp",
         representation = representation(
             def = "expression",
             covAlls = "list",
             hasGrad = "logical",
             label = "character",
             d = "integer",
             parN = "integer",
             parNames = "character", 
             inputNames = "character",      ## spatial var names length d
             topParN = "integer",
             topParNames = "character",
             topPar = "numeric",
             topParLower = "numeric",
             topParUpper = "numeric",
             parsedFormula = "list"
             ))
         ## contains = "covAll")

##============================================================================
##' Creator for the class "covComp" for Composite Covariance kernels.
##'
##' A covariance object is built using \code{formula} which involves
##' kernel objects inheriting from the class \code{"covAll"} and
##' possibly of other scalar numeric parameters called \emph{top}
##' parameters. The formula can be thought of as involving the
##' covariance matrices rather than the kernel objects, each kernel
##' object say \code{obj} being replaced by \code{covMat(obj, X)} for
##' some design matrix or data frame \code{X}. Indeed, the sum or the
##' product of two kernel objects lead to a covariance which is simply
##' the sum or product of the kernel covariances. The top parameters
##' are considered as parameters of the covariance structure, as well
##' as the parameters of the covariance objects used in the
##' formula. Their value at the creation time will be used and thus
##' will serve as initial value in estimation.
##'
##' @title Creator for the Class \code{"covComp"} for Composite Covariance
##' Kernels.
##'
##' @param formula A formula. See \bold{Examples}.
##'
##' @param where An environment where the covariance kernels objects
##' and top parameters will be looked for.
##'
##' @param topParLower A numeric vector of lower bounds for the "top"
##' parameters.
##'
##' @param topParUpper A numeric vector of upper bounds for the "top"
##' parameters.
##'
##' @param trace Integer level of verbosity.
##'
##' @param ... Not used yet. For passing other slot values.
##'
##' @return An object with S4 class \code{"covComp"}.
##'
##' @section Caution: The class definition and its creator are to
##' regarded as a DRAFT, many changes being necessary until a stable
##' implementation will be reached. The functions relating to this
##' class are not for final users of GP models, but rather to those
##' interested in the conception and specification in view of a future
##' release of the \bold{kergp} package.
##'
##' @examples
##' ## =========================================================================
##' ## build some kernels (with their inputNames) in the global environment
##' ## =========================================================================
##' 
##' myCovExp3 <- kMatern(d = 3, nu = "1/2")
##' inputNames(myCovExp3) <- c("x", "y", "z")
##'
##' myCovGauss2 <- kGauss(d = 2)
##' inputNames(myCovGauss2) <- c("temp1", "temp2")
##'
##' k <- kMatern(d = 1)
##' inputNames(k) <- "x"
##'
##' ell <- kMatern(d = 1)
##' inputNames(ell) <- "y"
##' 
##' tau2 <- 100
##' sigma2 <- 4
##' 
##' myCovComp <- covComp(formula = ~ tau2 * myCovGauss2() * myCovExp3() + sigma2 * k())
##' 
##' myCovComp1 <- covComp(formula = ~ myCovGauss2() * myCovExp3() + k())
##'
##' inputNames(myCovComp)
##' coef(myCovComp)
##' 
##' n <- 5
##' set.seed(1234)
##' X <- data.frame(x = runif(n), y = runif(n), z = runif(n),
##'                 temp1 = runif(n), temp2 = runif(n))
##' 
##' C <- covMat(myCovComp, X = X)
##'
##' Cg <- covMat(myCovComp, X = X, compGrad = TRUE)
##'
##' ## Simulation: purely formal example, not meaningful.
##' 
##' Y <- simulate(myCovComp, X = X, nsim = 100)
##' 
covComp <- function(formula,
                    where = .GlobalEnv,
                    ## topPar = NULL,
                    topParLower = NULL,
                    topParUpper = NULL,
                    trace = 0,
                    ...) {
        
    pf <- parseCovFormula(formula, where = where, trace = trace)

    hasGrad <- all(pf$hasGrad)
    
    covAlls <- list()
    
    for (kn in pf$kernNames) {
        covAlls[[kn]] <- get(kn, envir = where)
    }

    topParNames <- pf$parNamesList[[".top"]]
    topParN <- length(topParNames)
    
    if (topParN) {
        topPar <- rep(NA, topParN)
        names(topPar) <- topParNames
        
        ## retrieve top parameters from the environment 'where'
        for (pn in pf$parNamesList[[".top"]]) {
            res <- try(get(pn, envir = where, mode = "numeric"))
            if (!is(res, "try-error")) {
                topPar[[pn]] <- res
            }
        }
    } else {
        topPar <- numeric(0)
    }
     
    ## user-defined values ???
    ## topPar <- checkPar(topPar, parN = topParN, parNames = topParNames,
    ##                    default = 1.0)
    topParLower <- checkPar(topParLower, parN = topParN,
                            parNames = topParNames, default = -Inf)
    topParUpper <- checkPar(topParUpper, parN = topParN,
                            parNames = topParNames, default = Inf)    
    new("covComp",
        def = pf$expr,
        covAlls = covAlls,
        hasGrad = hasGrad,
        label = "",
        d = pf$d,      
        inputNames = pf$inputNames,      ## spatial var names length d
        parN = length(pf$parNames),
        parNames = pf$parNames,
        topParN = topParN,
        topParNames = topParNames,
        topPar = topPar,
        topParLower = topParLower,
        topParUpper = topParUpper,
        parsedFormula = pf
        )
    
}

setMethod("show", 
          signature = signature(object = "covComp"),
          definition = function(object){
              cat("'Composite' covariance kernel from expression\n")
              print(object@def)
              cat("\n")
              if (object@hasGrad) {
                  cat("o Analytical gradient is provided.\n\n")
              }
              cat("o Inputs required.\n")
              print(object@inputNames)
              cat("\n")
              
              cat("o Kernel (covAll objects).\n")
              print(names(object@covAlls))
              cat("\n")
              
              cat("o Parameters.\n")
              print(coef(object, as = "list"))
          })


setMethod("checkX",
          signature = signature(object = "covComp", X = "data.frame"),
          definition = function(object, X, strict = FALSE, ...){
              iN <- inputNames(object)
              if (strict) {
                  if (length(iN) != ncol(X) || !all(iN == colnames(X)) )
                      stop("colnames(X) and inpuNames(object) must be identical")
              }
              if (!all(inputNames(object) %in% colnames(X)))
                  stop("all the elements of inputNames(object) must be ",
                       " in colnames(X)")
              
              X <- X[ , iN, drop = FALSE]
          })

## **********************************************************************
## code for the covMat method
##
## **********************************************************************

.covMat_covComp <- function(object, X, Xnew = NULL,
                           compGrad = compGrad, checkNames = NULL, index = -1L) {    
    
    isXnew <- !is.null(Xnew)
    if (isXnew) compGrad <- FALSE
    
    pf <- object@parsedFormula
    
    for (kn in pf$kernNames) {
        val <-  covMat(object@covAlls[[kn]], X = X, Xnew = Xnew,
                       compGrad = compGrad)
        assign(x = kn, value = val)
        if (compGrad) {
            assign(x = paste(".", kn, sep = ""), value = attr(val, "gradient"))
        }
    }
    
    for (tpn in pf$topParNames) {
        val <- pf$topPar[tpn]
        assign(x = tpn, value = val)
    }
            
    ## C <- eval(parse(text = pf$simpleCall))
    C <- eval(pf$parsedSimpleCall)
    
    if (compGrad) {
        if (isXnew) {
            stop("'compGrad = TRUE' is not implemented for the case where ",
                 "'Xnew' is not NULL")
        }
        n <- nrow(X)
        parNames <- pf$parNames
        np <- length(parNames)
        dC <- array(NA,
                   dim = c(n, n, np),
                    dimnames = list(NULL, NULL, parNames))        
        for (pn in parNames) {
            ## dC[ , , pn] <- eval(parse(text = pf$gradExp[[pn]]))
            dC[ , , pn] <- eval(pf$parsedGradExp[[pn]])
        }
        attr(C, "gradient") <- dC
    }
    
    C
    
}

## **********************************************************************
## Note that with an object of class "covComp", the design 'X' can be
## a matrix or a data frame, because some of the inputs can be factor.
## **********************************************************************

setMethod("covMat",
          signature = signature(object = "covComp"), 
          definition = function(object, X, Xnew,
              compGrad = hasGrad(object), checkNames = NULL, index = -1L) {
              
              .covMat_covComp(object = object, X = X, Xnew = Xnew,
                             compGrad = compGrad,
                             checkNames = checkNames,
                             index = index)
              
          })



## *********************************************************************
## Note that with an object of class "covComp", the design 'X' can be
## a matrix or a data frame, because some of the inputs can be factor.
## *********************************************************************

.varVec_covComp <- function(object, X,
                           compGrad = FALSE, checkNames = NULL, index = -1L) {    

    ## cat("varVec_covComp\n")
    
    pf <- object@parsedFormula
    
    for (kn in pf$kernNames) {
        val <-  varVec(object@covAlls[[kn]], X = X,
                       compGrad = compGrad)
        assign(x = kn, value = val)
        if (compGrad) {
            assign(x = paste(".", kn, sep = ""), value = attr(val, "gradient"))
        }
    }
    
    for (tpn in pf$topParNames) {
        val <- pf$topPar[tpn]
        assign(x = tpn, value = val)
    }
            
    ## V <- eval(parse(text = pf$simpleCall))
    V <- eval(pf$parsedSimpleCall)
    names(V) <- rownames(X)
    
    if (compGrad) {
       
        n <- nrow(X)
        parNames <- pf$parNames
        np <- length(parNames)
        dV <- array(NA,
                   dim = c(n, np),
                   dimnames = list(NULL, parNames))
        for (pn in parNames) {
            ## dV[ , pn] <- eval(parse(text = pf$gradExp[[pn]]))
            dV[ , pn] <- eval(pf$parsedGradExp[[pn]])
        }
        
        attr(C, "gradient") <- dV
    }
    
    V
    
}

setMethod("varVec",
          signature = signature(object = "covComp"), 
          definition = function(object, X,
              compGrad = FALSE, checkNames = NULL, index = -1L) {
              
              .varVec_covComp(object = object, X = X,
                             compGrad = compGrad,
                             checkNames = checkNames,
                             index = index)
              
          })
              

## ===============================================================================
## `coef` and `coef<-` methods
## ===============================================================================

setMethod("coef",
          signature = signature(object = "covComp"), 
          definition = function(object, type = "all", as = "vector"){
              
              if (as == "vector") { 
                  res <- numeric(0)
                  for (i in seq_along(object@covAlls)) {
                      res <- c(res, coef(object@covAlls[[i]])) 
                  }
                  res <- c(res, object@topPar)
                  names(res) <- object@parNames
              } else if (as == "list") {
                  res <- list()
                  for (kn in names(object@covAlls)) {
                      res[[kn]] <- coef(object@covAlls[[kn]])
                  }
                  res[[".top"]] <- object@topPar
              }
              res
          })

setMethod("coef<-", 
          signature = signature(object = "covComp", value = "numeric"),
          definition = function(object, type = "all", as = "vector",
              ..., value){
              if (type != "all" || as != "vector") {
                  stop("at the time only implemented values: \"all\" for 'type' ",
                       "and \"vector\" for 'as'")
              }
              if (length(value) != object@parN) {
                  stop(sprintf("'value' must have length %d", object@parN))
              }
              used <- 0
              for (i in seq_along(object@covAlls)) {
                  npi <- object@covAlls[[i]]@parN
                  ind <- used + (1L:npi)
                  coef(object@covAlls[[i]]) <- value[ind]
                  used <- used + npi
              }
              if (object@topParN) {
                  ind <- (used + 1L):object@parN 
                  object@topPar <- value[ind]
              }
              object
              
      })

setMethod("coefLower",
          signature = signature(object = "covComp"), 
          definition = function(object, type = "all", as = "vector"){
              
              if (as == "vector") { 
                  res <- numeric(0)
                  for (i in seq_along(object@covAlls)) {
                      res <- c(res, coefLower(object@covAlls[[i]])) 
                  }
                  res <- c(res, object@topParLower)
                  names(res) <- object@parNames
              } else if (as == "list") {
                  res <- list()
                  for (kn in names(object@covAlls)) {
                      res[[kn]] <- coefLower(object@covAlls[[kn]])
                  }
                  res[[".top"]] <- object@topParLower
              }
              res
          })

setMethod("coefUpper",
          signature = signature(object = "covComp"), 
          definition = function(object, type = "all", as = "vector"){
              
              if (as == "vector") { 
                  res <- numeric(0)
                  for (i in seq_along(object@covAlls)) {
                      res <- c(res, coefUpper(object@covAlls[[i]])) 
                  }
                  res <- c(res, object@topParUpper)
                  names(res) <- object@parNames
              } else if (as == "list") {
                  res <- list()
                  for (kn in names(object@covAlls)) {
                      res[[kn]] <- coefUpper(object@covAlls[[kn]])
                  }
                  res[[".top"]] <- object@topParUpper
              }
              res
          })

## ===============================================================================
## to prevent unwanted results
## ===============================================================================

setMethod("inputNames<-", 
          signature = signature(object = "covComp", value = "character"),
          definition = function(object, ..., value){
              warning("replacement method `inputNames<-` not yet available ",
                      "for the class \"covComp\". Please use this method on the ",
                      "component covariance objects instead.")
              object
        })

##***********************************************************************
## scores method 
## 
## Note that the scores method can only be used when the weight matrix
## is known, which requires an evaluation of 'covMat'.
##
## To avoid that, the covariance matrix could be passed to the 'scores'
## method, whith a previous analysis of the kind of object (whether
## it accepts matrix or not)
##
##***********************************************************************

setMethod("scores",
          signature = "covComp", 
          definition = function(object, X, weights,  ...) {

              C <- covMat(object, X, checkNames = FALSE, compGrad = TRUE)
              dCov <- attr(C, "gradient")
              
              n <- nrow(X)
              if (any(dim(dCov) != c(n, n, object@parN))) {
                  stop("\"gradient\" attribute with wrong dimension")  
              }
              lt <- lower.tri(matrix(NA, nrow = n , ncol = n), diag = TRUE)
              agFun <- function(mat) sum(weights * mat[lt])
              scores <- apply(dCov, MARGIN = 3L, FUN = agFun)
              return(scores)
              
          })

##***********************************************************************
## Coercion into a list: useful to get the kernels after an estimation.
##***********************************************************************

setMethod(f = "as.list",
          signature = signature(x = "covComp"),
          definition = function(x) x@covAlls)
