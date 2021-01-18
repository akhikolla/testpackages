
setClass("covQual",   	
         representation(
             covLevels = "function",         ## correlation or covariance for all levels
                                             ## only one 'par' argument!
             covLevMat = "matrix",           ## added on 2017-10-11
             hasGrad = "logical",            ##
             acceptLowerSQRT = "logical",    ## if TRUE, we can pass 'lowerSQRT' to the fun
             label = "character",
             d = "integer",                  ## (spatial) dimension
             inputNames = "character",       ## spatial var names length d
             nlevels = "integer",            ##        
             levels = "list",                ## list of levels (labels)
             parLower = "numeric",           ## lower bound on pars
             parUpper = "numeric",           ## upper bound on pars
             par  = "numeric",               ## params values
             parN = "integer",               ## number of par
             kernParNames  = "character",    ## depending on kernel
             ordered = "logical",
             intAsChar = "logical"
         ),
         validity = function(object){
             if (length(object@kernParNames) != object@parN) {
                 stop("Incorrect number of parameter names")
             }
         })
         ## contains = "covAll")

## ***********************************************************************
## XXX DOES NOT WORK: reset class union. Is this really correct???
## This will not anymore a problem when 'catgp' will duely be
## integrated in 'kergp'.
## ***********************************************************************
## 
## setClassUnion("covAll", c("covTS", "covMan", "covQual"))

setGeneric("checkX",
           function(object, X, ...) standardGeneric("checkX")
           )

##  npar method for class "covQual".
setMethod("npar",
          signature = signature(object = "covQual"),
          definition = function(object,  ...){
            object@parN
          })

## setMethod("sd2",
##           signature = signature(object = )

##***********************************************************************
## CAUTION:  when 'type' is a vector and 'as' is "list" or "matrix"
## elements are returned in the order given by 'type'
## which might differ from the standard parameter order.
##
## o 'type' can be "all", or can be a character vector describing a
##          subset of the union U(kernParNaems, "var")
## 
## o 'as'   can be "vector", "list", or "matrix"
##
##***********************************************************************
setMethod("coef", 
          signature = signature(object = "covQual"), 
          definition = function(object){         
            res <- object@par
            names(res) <- object@kernParNames
            res
          })


##***********************************************************************
## Replacement method
##
## XXX check validity???
##
## NOT WRITTEN YET
##
##**********************************************************************

setMethod("coef<-", 
          signature = signature(object = "covQual", value = "numeric"),
          definition = function(object, ..., value){
            if (length(value) != object@parN) {
              stop(sprintf("'value' has length %d but must have length %d",
                           length(value), object@parN))
            }
            object@par[] <- value
            ## replace the slot
            object@covLevMat <- object@covLevels(value)
            colnames(object@covLevMat) <- rownames(object@covLevMat) <-
                object@levels[[1]]
            object
          })

##***********************************************************************
## Methods to get/set the parameter bounds?
## One could set bounds by group: range, shape etc.
##
##***********************************************************************
setMethod("coefLower", 
          signature = signature(object = "covQual"),
          definition = function(object){
              res <- object@parLower
              names(res) <- object@kernParNames
              res         
          })

setMethod("coefLower<-",
          signature = signature(object = "covQual"),
          definition = function(object, ..., value){
            if (length(value) != object@parN) {
              stop(sprintf("'value' must have length %d", object@parN))
            }
            object@parLower[] <- value
            object   
          })

setMethod("coefUpper", 
          signature = signature(object = "covQual"),
          definition = function(object){
              res <- object@parUpper
              names(res) <- object@kernParNames
              res            
          })

setMethod("coefUpper<-",
          signature = signature(object = "covQual"),
          definition = function(object, ..., value){
            if (length(value) != object@parN) {
              stop(sprintf("'value' must have length %d", object@parN))
            }
            object@parUpper[] <- value
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
          signature = "covQual", 
          definition = function(object, X, weights,  ...) {

              C <- covMat(object, X, checkNames = FALSE, compGrad = TRUE)
              dCov <- attr(C, "gradient")
              
              n <- nrow(X)
              if (any(dim(dCov) != c(n, n, npar(object)))) {
                  stop("\"gradient\" attribute with wrong dimension")  
              }
              lt <- lower.tri(matrix(NA, nrow = n , ncol = n), diag = TRUE)
              agFun <- function(mat) sum(weights * mat[lt])
              scores <- apply(dCov, MARGIN = 3L, FUN = agFun)
              return(scores)
              
              
          })


##***********************************************************************
## The 'show' method must show the kernel name and parameter structure.
## It should also provide information of the parameterisation of the
## structure itself (sharing of the parameters across inputs).
##
##***********************************************************************
setMethod("show", 
          signature = signature(object = "covQual"), 
          definition = function(object){
            cat("Qualitative covariance kernel\n")
            argNames <- names(formals(object@covLevels))
            cat(paste("o Description:"), object@label, "\n")
            cat(sprintf("o Dimension 'd' (nb of inputs): %d\n", object@d))
#             cat(paste("o Kernel depending on: \"", 
#                       argNames[1], "\", \"", argNames[2], "\"\n", sep=""))
            cat(paste("o Parameters: ",
                      paste(sprintf("\"%s\"", object@kernParNames),
                              collapse = ", "),
                      "\n",sep = ""))
            cat(sprintf("o Number of parameters: %d\n", object@parN))
            if (object@acceptLowerSQRT) {
                cat("o Accept Lower SQRT.\n")
            }
            if (object@hasGrad) {
                cat("o Analytical gradient is provided.\n")
            }
            cat("o Param. values: \n")
            co <- cbind(coef(object), coefLower(object), coefUpper(object))
            colnames(co) <- c("Value", "Lower", "Upper")
            print(co)
            if (object@d == 1L) {
                if (object@nlevels[1] <= 10) {
                    cat("o Covariance between levels: \n")
                    M <- covMat(object)
                    attr(M, "gradient") <- NULL
                    print(M)
                } else {
                    cat("use 'covMat' with only one argument to see the ",
                        "covariance matrix.\n")
                }
            }       
          })

## *************************************************************************
## covMat method
## *************************************************************************
setMethod("covMat",
          signature = signature(object = "covQual"), 
          definition = function(object, X, Xnew = NULL,
              compGrad = hasGrad(object), 
              checkNames = NULL,
              intAsChar = object@intAsChar,
              lowerSQRT = FALSE,
              ...) {
              
              isXnew <- !is.null(Xnew)
              if (isXnew) compGrad <- FALSE
              
              if (missing(X)) {
                  
                  ## this exit was added on 2017-10-11
 
                  # if (!compGrad && !lowerSQRT) {
                  if (!lowerSQRT) {
                      return(object@covLevMat)
                  }

                  if (object@d == 1L) {
                      ## CAUTION without levels = this would not
                      ## work as expected for ordinal factors.
                      f <- factor(seq_along(object@levels[[1]]),
                                  labels = object@levels[[1]])
                      X <- data.frame(f)
                      colnames(X) <- object@inputNames[1L]
                      rownames(X) <- object@levels[[1]]
                  } else {
                      stop("'X = NULL' is only allowed when 'object' ",
                           "has only one input")
                  }
              }
              
              ## X <- as.matrix(X)
              
              if (is.null(checkNames)) {
                  checkNames <- TRUE 
                  # if (object@d == 1L) {
                  #     if (ncol(X) == 1L) {
                  #         checkNames <- FALSE
                  #     }  
                  # }      
              }
              
              if (checkNames) X <- checkX(object, X = X, intAsChar = intAsChar)
              if (any(is.na(X))) stop("'X' must not contain NA elements")

              if (isXnew){
                  ## Xnew <- as.matrix(Xnew)
                  if (checkNames) Xnew <- checkX(object, X = Xnew, intAsChar = intAsChar)
                  if (ncol(X) != ncol(Xnew)) {
                      stop("'X' and 'Xnew' must have the same number of columns")
                  }
                  if (any(is.na(Xnew))) stop("'Xnew' must not contain NA elements") 
              } else {
                  Xnew <- X
              }
              
              compGrad <- as.integer(compGrad)
              par <- coef(object)  
              covLevs <- object@covLevels(object@par,
                                          lowerSQRT = lowerSQRT,
                                          compGrad = compGrad)
              if (object@d > 1L) {
                  stop("only one factor is allowed for now")
              }
              iX <- as.integer(X[ , 1L])
              if (!isXnew){
                  Cov <- covLevs[iX, iX, drop = FALSE]
              } else {
                  iXnew <- as.integer(Xnew[ , 1L])
                  if (compGrad) {
                      stop("Gradient computation not implemented when Xnew != NULL")
                  }
                  Cov <- covLevs[iX, iXnew, drop = FALSE]
              }
              rownames(Cov) <- rownames(X)
              colnames(Cov) <- rownames(Xnew)
              
              if (compGrad) {
                  ## do an apply
                  gradLevs <- attr(covLevs, "gradient")
                  ## print(gradLevs)
                  n <- nrow(X)
                  np <- npar(object)

                  ## this is best done with an 'apply' later, so not used
                  ##     grad <- array(NA, dim = c(n, n, np),
                  ##                   dimnames = list(rownames(Cov),
                  ##                       colnames(Cov), object@kernParNames))
                  ##     for (k in 1L:np) grad[ , , k] <- gradLevs[X, X, k] 
                                
                  grad <- apply(gradLevs, MARGIN = 3L, FUN = function(A) A[iX, iX])
                  dim(grad) <- c(n, n, np)
                  dimnames(grad) <- list(rownames(Cov),
                                         colnames(Cov),
                                         object@kernParNames)
                  attr(Cov, "gradient") <- grad
              }
                  
              return(Cov) 
              
          })

## *************************************************************************
## varVec method: compute the variance vector.
## *************************************************************************
setMethod("varVec",
          signature = signature(object = "covQual"), 
          definition = function(object, X,
              compGrad = FALSE,
              checkNames = NULL,
              intAsChar = object@intAsChar,
              ...) {
              
              if (missing(X)) {
                  if (object@d == 1L) {
                      f <- factor(object@levels[[1]])
                      X <- data.frame(f)
                      colnames(X) <- object@inputNames[1L]
                      rownames(X) <- object@levels[[1]]
                  } else {
                      stop("missing 'X' is only allowed when 'object' ",
                           "has only one input")
                  }
              }
              
              if (is.null(checkNames)) {
                  checkNames <- TRUE 
                  if (object@d == 1L) {
                      if (ncol(X) == 1L) checkNames <- FALSE
                  }      
              }
              
              if (checkNames) X <- checkX(object, X = X, intAsChar = intAsChar)
              if (any(is.na(X))) stop("'X' must not contain NA elements")
              
              compGrad <- as.integer(compGrad)
              par <- coef(object)  
              covLevs <- object@covLevels(object@par,
                                          lowerSQRT = FALSE,
                                          compGrad = compGrad)
              if (object@d > 1L) {
                  stop("only one factor is allowed for now")
              }

              iX <- as.integer(X[ , 1L])
              Cov <- covLevs[iX, iX, drop = FALSE]
              Var <- diag(Cov)
              names(Var) <- rownames(X)
              
              if (compGrad) {
                  ## do an apply
                  gradLevs <- attr(covLevs, "gradient")
                  ## print(gradLevs)
                  n <- nrow(X)
                  np <- npar(object)

                  ## this is best done with an 'apply' later, so not used
                  ##     grad <- array(NA, dim = c(n, np),
                  ##                   dimnames = list(rownames(X),
                  ##                                   object@kernParNames))
                  ##     for (k in 1L:np) grad[ , k] <- diag(gradLevs[iX, iX, k]) 
                                
                  grad <- apply(gradLevs, MARGIN = 3L,
                                FUN = function(A) diag(A[iX, iX]))
                  dim(grad) <- c(n, np)
                  dimnames(grad) <- list(rownames(X), object@kernParNames)
                  attr(Var, "gradient") <- grad
              }
                  
              return(Var) 
              
          })

## *************************************************************************
## simulate from a covariance structure.
## *************************************************************************
setMethod("simulate",
          signature = signature(object = "covQual"),
          definition = function(object,  nsim = 1,
              seed = NULL, X, mu = NULL, method = "mvrnorm",
              checkNames = TRUE, ...) {
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

              if (!missing(X)) {
                  if (iD <- anyDuplicated(X)) {
                      stop("duplicated value in row ", iD, ". The simulation ",
                           "of a stochastic process at duplicated indices is ",
                           "meaningless.")
                  }
                  C <- covMat(object, X = X, checkNames = checkNames)
              } else C <- covMat(object, checkNames = checkNames)
              
              n <- ncol(C)
              if (is.null(mu)) {
                  mu <- rep(0, n)
              } else {
                  mu <- rep(mu, length.out = n)
              }
              val <- t(mvrnorm(n = nsim, mu = mu, Sigma = C, ...))
              ## attr(val, "seed") <- RNGstate
              val
          })



## *************************************************************************
## plot a covariance structure.
## *************************************************************************

plot.covQual <- function(x, type = c("cov", "cor", "warping"), ...){
  type <- match.arg(type)
  if (type %in% c("cor", "cov")){
    C <- covMat(x)
    is.corr <- (type == "cor")
    if (type == "cor") C <- cov2cor(C)
    if (requireNamespace("corrplot", quietly = TRUE)){
      corrplot::corrplot(C, is.corr = is.corr, ...)
    } else {
      levelplot(C, ...)
      # image(C[nrow(C):1, nrow(C):1], xaxt = "n", yaxt = "n", ...)
      # axis(side = 3, at = seq(0, 1, length.out = x@nlevels), labels = x@levels[[1]])
      # axis(side = 2, at = seq(0, 1, length.out = x@nlevels), labels = rev(x@levels[[1]]))
    }
  }
  if (type == "warping"){
    if (class(x) != "covOrd"){
      warning("Type 'warping' only works for an ordinal factor")
      return(NULL)
    }
    parWarp <- coef(x)[1:x@parNwarp]
    zseq <- seq(0, 1, length.out = 200)
    yseq <- x@warpFun$fun(zseq, par = parWarp, L = x@nlevels)
    plot(zseq, yseq, type = "l", xaxt = "n", 
         xlab = inputNames(x), ylab = "", ...)
    levseq <- seq(0, 1, length.out = x@nlevels)
    axis(side = 1, at = levseq, labels = x@levels[[1]])
    points(levseq, x@warpFun$fun(levseq, par = parWarp, L = x@nlevels), pch = 19)
    if (x@parNk1 > 0) lines(c(0, 1), c(0, 1), lty = "dotted")  # in that case, F is normalized
  }
}

if(!isGeneric("plot")) {
  setGeneric(name = "plot",
             def = function(x, y, ...) standardGeneric("plot")
  )
}

setMethod("plot",
          signature(x = "covQual"), 
          function(x, y, type = c("cov", "cor", "warping"), ...){
            if (!missing(y)) warning("Argument y is ignored (not used in this plot method)")
            plot.covQual(x = x, type = type, ...)
          }
          )

## **********************************************************************
## Check that the QUALITATIVE design matrix X is compatible with the
## covariance object This is quite different form the continuous input
## case, because one has to check the content of 'X', and not only the
## input names.
## **********************************************************************

setMethod("checkX",
          signature = signature(object = "covQual", X = "matrix"),
          definition = function(object, X, strict = FALSE, ...){
              
              iN <- inputNames(object)
              
              if (strict) {
                  if (length(iN) != ncol(X) || !all(iN == colnames(X)) )
                      stop("colnames(X) and inpuNames(object) must be identical")
              }
              if ( !all(inputNames(object) %in% colnames(X)) )
                  stop("all the elements of inputNames(object) must be ",
                       "in colnames(X)")
              
              X <- X[ , iN, drop = FALSE]
              nX <- ncol(X)

              checkX(object = object, X = as.data.frame(X))

          }
          )
    

setMethod("checkX",
         signature = signature(object = "covQual", X = "data.frame"),
         definition = function(object, X, strict = FALSE, intAsChar = object@intAsChar, ...){
             iN <- inputNames(object)
             if (strict) {
                 if (length(iN) != ncol(X) || !all(iN == colnames(X)) )
                     stop("colnames(X) and inpuNames(object) must be identical")
             }
             if ( !all(inputNames(object) %in% colnames(X)) )
                 stop("all the elements of inputNames(object) must be",
                      " in colnames(X)")
             u <- X[ , iN]
             levs <- object@levels[[1]]

             if (is.integer(u) && intAsChar) {
                 u <- as.character(u)
             }
             
             if (is.integer(u)) {
                 if (any(u < 1) || any(u > object@nlevels)) {
                     stop("When 'inputNames(object)' corresponds to an ",
                          "integer column in 'X', this must contain values ",
                          "between 1 and the number of levels of 'object'")
                 }
                 if (object@ordered[1]) {  # works only for one factor
                   u <- ordered(levs[u], levels = levs)
                   warning("Automatic coercion used when input levels are integer") 
                   #         The first observations are: ")
                   # print(u[1:min(4, nlevels(object)[1]])
                 } else {
                   u <- factor(levs[u], levels = levs)
                 }
             } else if (is.character(u)) {
                 if (!all(u %in% levs)) {
                     stop("When 'inputNames(object)' corresponds to a ",
                          "character column in 'X', this must contain only ",
                          "values found in the levels of 'object'")
                 }
                 names(levs) <- levs
                 if (object@ordered[1]) {  # works only for one factor
                    u <- ordered(levs[u], levels = levs)
                 } else {
                    u <- factor(levs[u], levels = levs)
                 }
               
             }  else if (is.factor(u)) {

                 levsu <- levels(u)
                 if (!identical(levsu, levs)) {        
                     if (!setequal(levsu, levs)) {
                         stop("When 'inputNames(object)' corresponds to a ",
                              "factor column in 'X', the levels must be identical ",
                              "to the levels of 'object'")
                     } else {
                         warning("The levels of the factor input are not in the same ",
                                 "order as that of the 'covQual' object. Forcing the same ",
                                 "order")
                         u <- factor(u, levels = levs)
                     }
                 }
                 
                 if (object@ordered[1]) {  # works only for one factor
                   u <- as.ordered(u)
                 } 
                 
             } else  {
                 stop("'inputNames(object)' must correspond to a column in 'X' ",
                      "which is 'integer', 'character' or 'factor'")
             }
             df <- data.frame(u)
             colnames(df) <- iN
             df
           })
