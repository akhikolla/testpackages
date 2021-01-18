## ------------
## Union Class
## ------------

setClassUnion("covAll",
              c("covTS", "covMan", "covComp", "covQual", "covRadial"))


##***********************************************************************
## 'inputNames'
##
## XXX   This method should be removed for classes inheriting from
##       "covKernel", except if inputnames are stored under another
##       name.
## 
##***********************************************************************
setMethod("hasGrad",
          signature = signature(object = "covAll"),
          definition = function(object, ...){
            if ("hasGrad" %in% slotNames(object)) {
              object@hasGrad
            } else return(FALSE)
          })

setMethod("inputNames",
          signature = signature(object = "covAll"),
          definition = function(object, ...){
            object@inputNames
          })

setMethod("inputNames<-", 
          signature = signature(object = "covAll", value = "character"),
          definition = function(object, ..., value){
            if (length(value) != object@d) {
              stop(sprintf("'value' has length %d but must have length %d",
                           length(value), object@d))
            }
            object@inputNames[] <- value
            object
          })

##***********************************************************************
## Check that the design matrix X is compatible with the covariance
## object. 
##***********************************************************************

setMethod("checkX",
          signature = signature(object = "covAll", X = "matrix"),
          definition = function(object, X, strict = FALSE, ...){
            iN <- inputNames(object)
            if (strict) {
              if (length(iN) != ncol(X) || !all(iN == colnames(X)) )
                stop("colnames(X) and inpuNames(object) must be identical")
            }
            if ( !all(inputNames(object) %in% colnames(X)) )
                stop("all the elements of inputNames(object) must be",
                     " in colnames(X)")
            X[ , iN, drop = FALSE] 
          })

setMethod("checkX",
          signature = signature(object = "covAll", X = "data.frame"),
          definition = function(object, X, strict = FALSE, ...){
            iN <- inputNames(object)
            if (strict) {
              if (length(iN) != ncol(X) || !all(iN == colnames(X)) )
                stop("colnames(X) and inpuNames(object) must be identical")
            }
            if ( !all(inputNames(object) %in% colnames(X)) )
                stop("all the elements of inputNames(object) must be",
                     " in colnames(X)")
            X[ , iN, drop = FALSE] 
          })


##' Draw random values for the parameters of a covariance kernel
##' object.
##'
##' Draw random values for the parameters of a covariance kernel
##' object using the informations \code{coefLower} and
##' \code{coefUpper}.
##' 
##' @title Draw random values for the parameters of a covariance kernel.
##'
##' @param object A covariance kernel.
##'
##' @param nsim Number of drawings.
##'
##' @param seed Seed for the random generator.
##'
##' @param ... Other arguments for methods.
##'
##' @return A matrix with \code{nsim} rows and \code{npar(object)} columns.
##' 
##' @author ODY

setMethod("simulPar", 
          signature = signature(object = "covAll"),
          definition = function(object, nsim = 1L, seed = NULL){
            
            ## copied from simulate.lm
            if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
              runif(1)                     # initialize the RNG if necessary
            if(is.null(seed))
              RNGstate <- get(".Random.seed", envir = .GlobalEnv)
            else {
              R.seed <- get(".Random.seed", envir = .GlobalEnv)
              set.seed(seed)
              RNGstate <- structure(seed, kind = as.list(RNGkind()))
              on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
            }
            Par <- coef(object)
            L <- coefLower(object)
            U <- coefUpper(object)
            p <- length(Par)
            Par <- array(NA, dim = c(nsim, p), dimnames = list(NULL, names(Par)))
            for (k in  1:p) {
              if ( is.finite(L[k]) ) {
                if (is.finite(U[k])) {
                  Par[ , k] <- runif(n = nsim, min = L[k], max = U[k])
                } else {
                  Par[ , k] <- L[k] + rexp(n = nsim, rate = 1)
                } 
              } else{
                if ( is.finite(U[k]) ) {
                  Par[ , k] <- U[k] - rexp(n = nsim, rate = 1)
                } else {
                  Par[ , k] <- rcauchy(n = nsim)
                }   
              }
            }
            Par 
          }
          )



##***********************************************************************
## The 'simulate' method 
##
##' show method for class "covAll"
##' @aliases simulate,covAll-method
##'
##' @param object XXX
##' @docType methods
##' @rdname covAll-methods
##'
##***********************************************************************
## if (!isGeneric("simulate")) {
##   setGeneric(name = "simulate",
##              def = function(object, nsim = 1, seed = NULL, ...) standardGeneric("simulate")
##              )
## }

setMethod("simulate",
          signature = signature(object = "covAll"),
          definition = function(object,  nsim = 1L,
              seed = NULL, X, mu = NULL, method = "mvrnorm",
              checkNames = TRUE,
              ...) {
              require(MASS)
              ## X <- as.matrix(X)
              ## lines copied from 'simulate.lm'
              if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
                  runif(1)
              if (is.null(seed)) 
                  RNGstate <- get(".Random.seed", envir = .GlobalEnv)
              else {
                  R.seed <- get(".Random.seed", envir = .GlobalEnv)
                  set.seed(seed)
                  RNGstate <- structure(seed, kind = as.list(RNGkind()))
                  on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
              }
              C <- covMat(object, X = X, checkNames = checkNames, compGrad = FALSE)
              n <- ncol(C)
              if (is.null(mu)) {
                  mu <- rep(0, n)
              } else {
                  mu <- rep(mu, length.out = n)
              }

              if (nsim == 1L) {
                  val <- as.matrix(mvrnorm(n = nsim, mu = mu, Sigma = C, ...))
              } else {
                  val <- t(as.matrix(mvrnorm(n = nsim, mu = mu, Sigma = C, ...)))
              }
              ## attr(val, "seed") <- RNGstate
              val
          })


## ************************************************************************
## the 'varVec' method computes (only) the diagonal of the covMat with
## the same inputs, This will be used for classes of covariance
## kernels such as kernels for qualitative variable which are neither
## of class "covMan" nor of class "covTS"
##
## CAUTION: this was commented here because the 'varVec' method can not
## call the 'covMat' method at this level since the later does not exist.
##
## ************************************************************************
## setMethod("varVec",
##           signature = signature(object = "covAll", X = "matrix"),
##           definition = function(object,  X,
##               compGrad = FALSE,
##               checkNames = NULL,
##               ...) {
              
##               if (checkNames) X <- checkX(object, X = X)

##               if (compGrad) {
##                   warning("'compGrad = TRUE' not yet allowed here. Coming soon.")
##               }
                  
##               ## checkNames must be turned to FALSE for the sake of speed
##               apply(X, MARGIN = 1,
##                     FUN = function(x, ...) covMat(object, matrix(x, nrow = 1L), ...),
##                     compGrad = FALSE, checkNames = checkNames)
              

##           })

## setMethod("varVec",
##           signature = signature(object = "covAll", X = "data.frame"),
##           definition = function(object,  X,
##               compGrad = FALSE,
##               checkNames = NULL,
##               ...) {
              
##               if (!is.null(checkNames) && checkNames) X <- checkX(object, X = X)

##               if (compGrad) {
##                   warning("'compGrad = TRUE' not yet allowed here. Coming soon.")
##               }
                  
##               ## checkNames must be turned to FALSE for the sake of speed
##               apply(X, MARGIN = 1,
##                     FUN = function(x, ...) covMat(object, matrix(x, nrow = 1L), ...),
##                     compGrad = FALSE, checkNames = checkNames)
              

##           })
