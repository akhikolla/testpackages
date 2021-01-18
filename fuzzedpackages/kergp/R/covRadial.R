setClass("covRadial",   	
         representation(
             k1Fun1 = "function",
             hasGrad = "logical",
             cov = "integer",               ## 0 : corr, 1 : homo
             iso = "integer",               ## iso = 0
             label = "character",
             d = "integer",                 ## (spatial) dimension
             inputNames = "optCharacter",   ## spatial var names length d
             parLower = "numeric",          ## lower bound on par's
             parUpper = "numeric",          ## upper bound on par's
             par  = "numeric",              ## params values
             parN1 = "integer",             ## nbr of par in kern1fun, usually 0
             parN = "integer",              ## number of par
             kern1ParNames = "character", 
             kernParNames  = "character"    ## depending on kernel
         ),
         validity = function(object){
           if (length(object@kernParNames) != object@parN) {
               stop("Incorrect number of parameter names")
           }
           
         }
         )

## *****************************************************************************
##' Creator for the class \code{"covRadial"}.
##'
##' When \code{k1Fun1} has more than one formal argument, the
##' arguments with position \code{> 1} are assumed to be parameters of
##' the model. Examples are functions with formals \code{function(x,
##' shape = 1.0)} or \code{function(x, alpha = 2.0, beta = 3.0)},
##' corresponding to vector of parameter names \code{c("shape")} and
##' \code{c("alpha", "beta")}.
##'
##' @title Creator for the Class \code{"covRadial"}
##' 
##' @param k1Fun1 A function of a scalar numeric variable. Not that
##' using a one-dimensional kernel here \emph{does not} warrant that
##' a positive semi-definite kernel results for any diemnsion \eqn{d}.
##'
##' @param cov Integer. The value \code{0L} corresponds to a
##' correlation kernel while \code{1L} is for a covariance kernel.
##' 
##' @param iso Integer. The value \code{1L} coresponds to an isotropic
##' covariance, with all the inputs sharing the same range value.
##'
##' @param hasGrad Integer or logical. Tells if the velue returned by
##' the function \code{k1Fun1} has an attribute named \code{"der"}
##' giving the derivative(s).
##' 
##' @param inputs Character. Names of the inputs.
##'
##' @param d Integer. Number of inputs. 
##' 
##' @param parNames. Names of the parameters. By default, ranges
##' are prefixed \code{"theta_"} in the non-iso case and the range is
##' names \code{"theta"} in the iso case.
##' 
##' @param par Numeric values for the parameters. Can be \code{NA}.
##'
##' @param parLower Numeric values for the lower bounds on the
##' parameters. Can be \code{-Inf}.
##'
##' @param parUpper Numeric values for the uppper bounds on the
##' parameters. Can be \code{Inf}.
##'
##' @param ... Other arguments passed to the method \code{new}.
##'
##' @return An object with class \code{"covRadial"}.
##' 
##' 
covRadial <- function(k1Fun1 = k1Fun1Gauss,
                      cov = c("corr", "homo"),
                      iso = 0,
                      hasGrad = TRUE,
                      inputs = NULL,
                      d = NULL,
                      parNames,
                      par = NULL,
                      parLower = NULL,
                      parUpper = NULL,
                      label = "Radial kernel",
                      ...) {

    L <- list(...)
    
    cov <- match.arg(cov)
    cov <- switch(cov, corr = 0, homo = 1)
    
    if (is.null(d)) {
        if (is.null(inputs)) {
            stop ("at least one of the formals 'd' and 'inputs' must be given")
        } else {
            d <- length(inputs)
        }
    } else {
        if (is.null(inputs)) {
            inputs <-  paste("x", 1:d, sep = "")
        } else {
            if (length(inputs) != d) {
                stop("'inputs' must be of length ", d)
            }
        }
    }
    
    k1Fun1 <- match.fun(k1Fun1)

    ## =========================================================================
    ## If 'k1Fun1' has more than one formal argument, the
    ## arguments having position > 1 become parameters
    ## for the kernl, with their name unchanged
    ## ! Caution: possible conflicts here for some names
    ## like "theta" or "sigma2".
    ## =========================================================================
    f <- as.list(formals(k1Fun1))
    f[[1]] <- NULL
    parN1 <- length(f)
    
    if (parN1) {

        if (parN1 > 1) {
            stop("For now, only one \"shape\" parameter is allowed!")
        }
        
        kern1ParNames <- names(f)
        
        if (length(grep("theta", kern1ParNames)) ||
            length(grep("sigma", kern1ParNames))) {
            stop("The names of the formals arguments of 'k1Fun1' ",
                 "are in conflict with default names  \"theta\" and ",
                 "\"sigma\". Please rename.")
        }
        ## Set the shape parameters to their default value
        par1 <- sapply(f, function(x) if (x == "") return(NA) else return(x))        
    } else {
        kern1ParNames <- character(0)
    }

    parN <- parN1
    
    if (iso) {
        kernParNames <- c(kern1ParNames, "theta")
        parN <- parN + 1L
    } else {
        kernParNames <- c(kern1ParNames, paste0("theta_", 1:d))
        parN <- parN + d
    }

    if (cov) {
        kernParNames <- c(kernParNames, "sigma2")
        parN <- parN + 1L
    }
    
    if (is.null(par)) {
        par <- as.numeric(rep(1.0, parN))
        if (parN1) par[1:parN1] <- par1 
    }
    
    if (is.null(parLower)) {
        parLower <- as.numeric(rep(0, parN))
        if (parN1) parLower[1:parN1] <- -Inf 
    }
    if (is.null(parUpper)) {
        parUpper <- as.numeric(rep(Inf, parN))
    }
    
    if (missing(d) & missing(inputs)) {
        stop("at least one of 'd' or 'inputs' must be provided")
    }
    if (length(inputs) != d) {
        stop("'d' must be equal to 'length(inputs)'")
    }

    ## XXX check that the function is zero. The extra
    ## arguments should have a default value
    if (k1Fun1(0) != 1.0) {
        stop("the function given in 'k1Fun1' must be such that ",
             "k1Fun1(0) = 1.0")
    }
    
    ## XXX check that the derivative is given 
    if (hasGrad) {
        
    }
    
    new("covRadial", 
        k1Fun1 = k1Fun1,
        hasGrad = as.logical(hasGrad),
        cov = as.integer(cov),
        iso = as.integer(iso),
        label = as.character(label),
        d = as.integer(d),
        inputNames = as.character(inputs),
        par = as.numeric(par),
        parLower = as.numeric(parLower),
        parUpper = as.numeric(parUpper),
        parN1 = as.integer(parN1),
        parN = as.integer(parN),
        kern1ParNames = as.character(kern1ParNames),
        kernParNames = as.character(kernParNames),
        ...)
    
} ## XXX TODO : check that the kernel has 2 (eg : brownian) or 3
## arguments (parameterised kernel).


##*****************************************************************************
## 'covMat' method
##
## When 'Xnew' is not given, i.e. in the symmetric case, efforts are
## made to exploit the symmetry. In particular, the univariate kernel
## function (which can be very costly in some cases) is invoked only
## (n - 1) * n / 2 times.
##
## See the document 'kergp Computing Details' to understand the
## computation of the gradient or of the derivatives.
## 
##******************************************************************************
setMethod("covMat",
          signature = "covRadial", 
          definition = function(object, X, Xnew,
              compGrad = hasGrad(object),
              deriv = 0,
              checkNames = NULL, ...) {
              
              isXNew <- !is.null(Xnew)
              if (isXNew) compGrad <- FALSE
              
              ## ==============================================================
              ## Specific checks related to 'deriv'.
              ## ==============================================================
              
              if (deriv > 0) {
                  if (compGrad) {
                      stop("For now, 'deriv' can be non-zero only when ",
                           "'compGrad' is FALSE")
                  }
                  if (!isXNew) {
                      stop("When 'Xnew' is not given, 'deriv' can only be 0") 
                  }
                  if (deriv > 1) {
                      stop("Second-order derivatives are not yet available")
                  }
              }
              
              ## ==============================================================
              ## The kernel has only continuous (numeric) inputs...
              ## ==============================================================

              ## Moved after checkNames, or the matrix may be of mode
              ## "character!!!
              ## X <- as.matrix(X)
              
              if (is.null(checkNames)) {
                  checkNames <- TRUE 
                  if (object@d == 1L) {
                      if (ncol(X) == 1L) checkNames <- FALSE
                  }      
              }
              
              ## ==============================================================
              ## check names and content of 'X', and of 'Xnew' if
              ## provided.
              ## ==============================================================
      
              if (checkNames) X <- checkX(object, X = X)
              if (any(is.na(X))) stop("'X' must not contain NA elements")
              X <- as.matrix(X)
              
              if (isXNew){
                  ## XXX unclear for now if 'Xnew' being a vector is
                  ## accepted or not
                  ##
                  ## Xnew <- as.matrix(Xnew)
                  if (checkNames) Xnew <- checkX(object, X = Xnew)
                  Xnew <- as.matrix(Xnew)
                  if (ncol(X) != ncol(Xnew)) {
                      stop("'X' and 'Xnew' must have the same number of columns")
                  }
                  if (any(is.na(Xnew))) {
                      stop("'Xnew' must not contain NA elements")
                  }
                  nNew <- nrow(Xnew)
              } 

              ## ===============================================================
              ## Some parameters. The vector 'theta' of ranges is
              ## recycled to have length 'd' in the 'iso' case
              ## ===============================================================
              
              d <- object@d
              n <- nrow(X)
              compGrad <- as.integer(compGrad)
              parN1 <- object@parN1
              parN <- object@parN
              par <- coef(object)
              nTheta <- ifelse(object@iso, 1L, d)
              theta <- par[parN1 + (1L:nTheta)]
              theta <- rep(theta, length.out = d)
              if (parN1) psi <- par[1:parN1]
                  
              if (object@cov) sigma2 <- par[parN1 + nTheta + 1L]
              else sigma2 <- 1.0
      
              if (compGrad) {

                  if (!object@hasGrad) {
                      stop("'object' does not allow the computation of the",
                           " gradient")
                  }
                  if (isXNew) {
                      stop("the computation of the gradient is not allowed",
                           " for now when 'Xnew' is given (only the symmetric",
                           " case is allowed)")
                  }
                  
                  sI <- symIndices(n)

                  ## ==========================================================
                  ## We need to store the quantities H2 used in the
                  ## gradient Note that using apply is MUCH SLOWER
                  ## than using the 'H2Tot' matrix(!)
                  ## ==========================================================

                  m <- (n - 1L) * n / 2
                  H2 <- array(NA, dim = c(m, d))
                  H2Tot <- array(0, dim = c(m, 1))
                  
                  for (ell in 1:d) {
                      H2[ , ell] <- ((X[sI$i, ell] - X[sI$j, ell])
                                     / theta[ell])^2
                      H2Tot <- H2Tot + H2[ , ell]
                  }
                  
                  ## R <- sqrt(apply(H2, MARGIN = 1, sum))
                  R <- sqrt(H2Tot)

                  ## ==========================================================
                  ## !Care for possible extra arguments passed as a vector
                  ## 'psi'
                  ## ==========================================================
                  
                  if (parN1) {
                      Kappa <- do.call(object@k1Fun1, list(R, psi))
                  } else {
                      Kappa <- do.call(object@k1Fun1, list(R))
                  }
                 
                  dKappa <- attr(Kappa, "der")[ , 1L]
                  
                  ## ==========================================================
                  ## We now make
                  ##
                  ## o 'KappaSym' to be a n * n symmetric matrix and
                  ## 
                  ## o 'GradSym' an array with dim c(n, n, partN) and
                  ## with symmetric n * n slices GradSym[ , , k].
                  ##
                  ## Note that the temporary array 'Grad' is needed
                  ## because it is not possible to use a two-indice
                  ## style of indexation for a 3-dimensional
                  ## array. 'Grad' keeps a zero diagonal all along the
                  ## loop.
                  ## ==========================================================
                  
                  KappaSym <- array(1.0, dim = c(n, n))
                  KappaSym[sI$kL] <- KappaSym[sI$kU] <- Kappa
                  
                  GradSym <- array(0.0, dim = c(n, n, parN))
                  dimnames(GradSym) <- list(NULL, NULL, object@kernParNames)
                  Grad <- array(0.0, dim = c(n, n))
                  
                  ## 'shape' parameters of the kernel
                  if (parN1 > 0L) {
                      for (k in 1L:parN1) {
                          Grad[sI$kL] <- Grad[sI$kU] <- sigma2 *
                              attr(Kappa, "gradient")[ , k]
                          GradSym[ , , k] <- Grad
                      }
                  }

                  ## range parameter(s) 
                  if (!object@iso) {
                      for (ell in 1L:d) {
                          Grad[sI$kL] <- Grad[sI$kU] <- -sigma2 * H2[ , ell] *
                              dKappa / theta[ell] / R
                          GradSym[ , , parN1 + ell] <- Grad
                      }
                  } else {
                      Grad[sI$kL] <- Grad[sI$kU] <- -sigma2 * R * dKappa /
                          theta[1L]
                      GradSym[ , , parN1 + 1L] <- Grad
                  }
                  
                  ## variance parameter Now the diagonal of the slice is 1.0
                  if (object@cov) {
                      GradSym[ , , parN]  <- KappaSym
                      KappaSym <- sigma2 * KappaSym
                  }
                  
                  attr(KappaSym, "gradient") <- GradSym
                  return(KappaSym)
                  
              } else {

                  ## ==========================================================
                  ## General (unsymmetric) case.
                  ## ==========================================================
                  
                  if (isXNew) {

                      if (deriv == 0) {

                          ## ==================================================
                          ## When the gradient is not required, we do
                          ## not need to store the matrices H2.
                          ## ==================================================
                          
                          H2 <- array(0.0, dim = c(n, nNew))
                      
                          for (ell in 1:d) {
                              H2 <- H2 +
                                  (outer(X = X[ , ell], Y = Xnew[ , ell], "-") /
                                       theta[ell])^2
                          }
                          
                          R <- sqrt(H2)
                          
                          ## ==================================================
                          ## !Care for possible extra arguments passed
                          ## as a vector 'psi'
                          ## ==================================================
                          
                          if (parN1) {
                              Kappa <- do.call(object@k1Fun1, list(R, psi))
                          } else {
                              Kappa <- do.call(object@k1Fun1, list(R))
                          }
                          
                          Kappa <- array(Kappa, dim = c(n, nNew))
                          if (object@cov) Kappa <- sigma2 * Kappa

                          return(Kappa)

                          
                      } else {
                          
                          ## ===================================================
                          ## We need to store the quantities H2 used
                          ## in the gradient
                          ## ===================================================
                          
                          H <- array(0.0, dim = c(n, nNew, d))
                          H2Tot <- array(0.0, dim = c(n, nNew))
                          
                          for (ell in 1:d) {
                              H[ , , ell] <-
                                  (outer(X[ , ell], Xnew[ , ell], "-") /
                                       theta[ell])
                              H2Tot <- H2Tot + H[ , , ell]^2
                          }
                  
                          ## R <- sqrt(apply(H, MARGIN = c(1L, 2L),
                          ##                 function(x) sum(x * x)))
                          
                          R <- sqrt(H2Tot)

                          ## ==================================================
                          ## !Care for possible extra arguments passed
                          ## as a vector 'psi'
                          ## ==================================================
                          
                          if (parN1) {
                              Kappa <- do.call(object@k1Fun1, list(R, psi))
                          } else {
                              Kappa <- do.call(object@k1Fun1, list(R))
                          }

                          dim(Kappa) <- c(n, nNew)
                          dKappa <- array(attr(Kappa, "der")[ , 1L],
                                          dim = c(n, nNew))
                          if (object@cov) {
                              Kappa <- sigma2 * Kappa
                          }
                          
                          ## ==================================================
                          ## ! It can be the case that 'R' contains
                          ## zeros, because a row of 'Xnew' can be
                          ## identical to a row of 'X'. 
                          ## ==================================================
                          
                          der <- array(0.0, dim = c(n, nNew, d))
                          dimnames(der) <- list(NULL, NULL, colnames(X))
                          eps <- sqrt(.Machine$double.eps)
                          
                          for (ell in 1L:d) {
                              prov <- H[ , , ell] / R
                              prov[R < eps] <- 0.0
                              der[ , , ell] <- - sigma2 * prov * dKappa /
                                  theta[ell]
                          }

                          ## note that 'kappa' already has a 'der' attribute! 
                          attr(Kappa, "der") <- der
                          return(Kappa)
                          
                      }
                          
                  } else {

                      sI <- symIndices(n)

                      ## ======================================================
                      ## We no longer need to store the quantities
                      ## 'H2' computed here as vectors with length
                      ## 'm'.
                      ## ======================================================
                      
                      m <- (n - 1L) * n / 2
                      H2 <- rep(0.0, m)
                      
                      for (ell in 1:d) {
                          H2 <- H2 + ((X[sI$i, ell] - X[sI$j, ell]) /
                                          theta[ell])^2
                      }
                      
                      R <- sqrt(H2)
                      
                      ## =======================================================
                      ## !Care for possible extra arguments passed as a vector
                      ## 'psi'
                      ## =======================================================
                      
                      if (parN1) {
                          Kappa <- do.call(object@k1Fun1, list(R, psi))
                      } else {
                          Kappa <- do.call(object@k1Fun1, list(R))
                      }
                                       
                      ## =======================================================
                      ## Reshape to a symmetric matrix. We assume that
                      ## the kernel function is a correlation function
                      ## so the diagonal of the returned matrix
                      ## contains ones.
                      ## =======================================================
                      
                      ## kf <- object@k1Fun1(0.0)
                      if (object@cov) Kappa <- sigma2 * Kappa
                      KappaSym <- array(sigma2, dim = c(n, n))
                      KappaSym[sI$kL] <- KappaSym[sI$kU] <- Kappa                    

                      return(KappaSym)
                      
                  }

              }
              
          })



## *************************************************************************
## 'varVec' method: compute the variance vector.
## *************************************************************************
setMethod("varVec",
          signature = "covRadial", 
          definition = function(object, X, compGrad = FALSE, 
                                checkNames = NULL, index = -1L, ...) {
              
              if (is.null(checkNames)) {
                  checkNames <- TRUE 
                  if (object@d == 1L) {
                      if (ncol(X) == 1L) checkNames <- FALSE
                  }      
              }
              
              if (checkNames) X <- checkX(object, X = X)
              if (any(is.na(X))) stop("'X' must not contain NA elements")
              
              if (compGrad) {
                  warning("'compGrad' is not yet implemented. Ignored.")
              }
              
              if (object@cov) {
                  Var <- rep(unname(coef(object)[object@parN]), nrow(X))
              } else {
                  Var <- rep(1.0, nrow(X))
              }
             
              
              return(Var) 
              
          })

## *************************************************************************
##' npar method for class "coRadial".
##'
##' npar method for the "covRadial" class
##'
##' @param object An object with class "covRadial"
##' 
##' @return The number of free parmaeters in a `covRadial`covariance.
##'
##' @docType methods
##'
##' @rdname covRadial-methods
##'
setMethod("npar",
          signature = signature(object = "covRadial"),
          definition = function(object,  ...){
            object@parN
          })


setMethod("coef", 
          signature = signature(object = "covRadial"), 
          definition = function(object){         
            res <- object@par
            names(res) <- object@kernParNames
            res
          })

setMethod("coef<-", 
          signature = signature(object = "covRadial", value = "numeric"),
          definition = function(object, ..., value){
            if (length(value) != object@parN) {
              stop(sprintf("'value' has length %d but must have length %d",
                           length(value), object@parN))
            }
            object@par[] <- value
            object
          })

##***********************************************************************
##***********************************************************************
setMethod("coefLower", 
          signature = signature(object = "covRadial"),
          definition = function(object){
            object@parLower            
          })

setMethod("coefLower<-",
          signature = signature(object = "covRadial"),
          definition = function(object, ..., value){
            if (length(value) != object@parN) {
              stop(sprintf("'value' must have length %d", object@parN))
            }
            object@parLower[] <- value
            object   
          })

setMethod("coefUpper", 
          signature = signature(object = "covRadial"),
          definition = function(object){
            object@parUpper            
          })

setMethod("coefUpper<-",
          signature = signature(object = "covRadial"),
          definition = function(object, ..., value){
            if (length(value) != object@parN) {
              stop(sprintf("'value' must have length %d", object@parN))
            }
            object@parUpper[] <- value
            object   
          })

##***********************************************************************
## coercion method to cleanly extract the kernel slot
## XXX to be written
##***********************************************************************

## setAs("covRadial", "function", function(from) from@kernel)

##***********************************************************************
## scores method 
## 
## Note that the scores method can only be used when the weight matrix
## is known.
##***********************************************************************
setMethod("scores",
          signature = "covRadial", 
          definition = function(object, X, weights, ...) {
              
              X <- as.matrix(X)
              n <- nrow(X)
              d <- ncol(X)
              if (any(is.na(X))) stop("'X' must not contain NA elements") 
              Cov <- covMat(object = object, X = X, compGrad = TRUE)
              
              dCov <- attr(Cov, "gradient")
              if (any(dim(dCov) != c(n, n, object@parN))) {
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
##' show method for class "covRadial"
##' @aliases show,covRadial-method
##'
##' @param object XXX
##' 
##' @docType methods
##'
##' @rdname covRadial-methods
##'
##***********************************************************************
setMethod("show", 
          signature = signature(object = "covRadial"), 
          definition = function(object){
            cat("Radial covariance kernel\n")
            cat(paste("o Description:"), object@label, "\n")
            cat(sprintf("o Dimension 'd' (nb of inputs): %d\n", object@d))
#             cat(paste("o Kernel depending on: \"", 
#                       argNames[1], "\", \"", argNames[2], "\"\n", sep=""))
            cat(paste("o Parameters: ",
                      paste(sprintf("\"%s\"", object@kernParNames),
                              collapse = ", "),
                      "\n",sep = ""))
            cat(sprintf("o Number of parameters: %d\n", object@parN))
            if (object@hasGrad) {
                cat("o Analytical gradient is provided.\n")
            }
            cat("o Param. values: \n")
            co <- cbind(coef(object), coefLower(object), coefUpper(object))
            colnames(co) <- c("Value", "Lower", "Upper")
            print(co)
          })



