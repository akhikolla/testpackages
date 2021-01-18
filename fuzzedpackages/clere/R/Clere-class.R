######################################################################
############################ Class Clere #############################
############################## Creation ##############################
######################################################################

#' \code{\linkS4class{Clere}} class
#' 
#' This class contains all the input parameters to run CLERE.
#' 
#' \describe{
#'   \item{y}{[numeric]: The vector of observed responses.}
#'   \item{x}{[matrix]: The matrix of predictors.} 
#'   \item{n}{[integer]: The sample size or the number of rows in matrix x.} 
#'   \item{p}{[integer]: The number of variables of the number of columns in matrix x.} 
#'   \item{g}{[integer]: The number or the maximum number of groups considered. Maximum number of groups stands when model selection is required.} 
#'   \item{nItMC}{[numeric]: Number of Gibbs iterations to generate the partitions.} 
#'   \item{nItEM}{[numeric]: Number of SEM/MCEM iterations.} 
#'   \item{nBurn}{[numeric]: Number of SEM iterations discarded before calculating the MLE which is averaged over SEM draws.}
#'   \item{dp}{[numeric]: Number of iterations between sampled partitions when calculating the likelihood at the end of the run.} 
#'   \item{nsamp}{[numeric]: Number of sampled partitions for calculating the likelihood at the end of the run.} 
#'   \item{sparse}{[logical]: Should a \code{0} class be imposed to the model?} 
#'   \item{analysis}{[character]: Which analysis is to be performed. Values are \code{"fit"}, \code{"bic"}, \code{"aic"} and \code{"icl"}.} 
#'   \item{algorithm}{[character]: The algorithmto be chosen to fit the model. Either the SEM-Gibbs algorithm or the MCEM algorithm. The most efficient algorithm being the SEM-Gibbs approach. MCEM is not available for binary response.} 
#'   \item{initialized}{[logical]: Is set to TRUE when an initial partition and an initial vector of parameters is given by the user.} 
#'   \item{maxit}{[numeric]: An EM algorithm is used inside the SEM to maximize the complete log-likelihood \code{p(y,Z|theta)}. \code{maxit} stands as the maximum number of EM iterations for the internal EM.} 
#'   \item{tol}{[numeric]: Maximum increased in complete log-likelihood for the internal EM (stopping criterion).} 
#'   \item{seed}{[integer]: An integer given as a seed for random number generation. If set to \code{NULL}, then a random seed is generated between \code{1} and \code{1000}.}
#'   \item{b}{[numeric]: Vector of parameter b. Its size equals the number of group(s).} 
#'   \item{pi}{[numeric]: Vector of parameter pi. Its size equals the number of group(s).} 
#'   \item{sigma2}{[numeric]: Parameter sigma^2.}
#'   \item{gamma2}{[numeric]: Parameter gamma^2.} itemintercept[numeric]: Parameter beta_0 (intercept).
#'   \item{likelihood}{[numeric]: Approximated log-likelihood.} 
#'   \item{entropy}{[numeric]: Approximated entropy.}
#'   \item{P}{[matrix]: A \code{p x g} matrix of posterior probability of membership to the groups. \code{P = E[Z|theta]}.} 
#'   \item{theta}{[matrix]: A \code{nItEM x (2g+4)} matrix containing values of the model parameters and complete data likelihood at each iteration of the SEM/MCEM algorithm}
#'   \item{Bw}{[matrix]: A \code{p x nsamp} matrix which columns are samples from the posterior distribution of Beta (regression coefficients) given the data and the maximum likelihood estimates.} 
#'   \item{Zw}{[matrix]: A \code{p x nsamp} matrix which columns are samples from the posterior distribution of Z (groups membership indicators) given the data and the maximum likelihood estimates.} 
#'   \item{theta0}{[numeric]: A \code{2g+3} length vector containing initial guess of the model parameters. See example for function \code{\link{fitClere}.}} 
#'   \item{Z0}{[numeric]: A \code{p x 1} vector of integers taking values between 1 and \code{p} (number of variables).}
#' }
#' 
#' @name Clere-class
#' @aliases Clere-class [,Clere-method [<-,Clere-method show,Clere-method sClere-class show,sClere-method [,sClere-method [<-,sClere-method
#' @docType class
#' 
#' @section Methods: 
#' \describe{ 
#'   \item{\code{object["slotName"]}:}{Get the value of the field \code{slotName}.} 
#'   \item{\code{object["slotName"]<-value}:}{Set \code{value} to the field \code{slotName}.} 
#'   \item{\code{plot(x, ...)}:}{Graphical summary for MCEM/SEM-Gibbs estimation.} 
#'   \item{\code{clusters(object, threshold = NULL, ...)}:}{Returns the estimated clustering of variables.}
#'   \item{\code{predict(object, newx, ...)}:}{Returns prediction using a fitted model and a new matrix of design.} 
#'   \item{\code{summary(object, ...)}:}{summarizes the output of function \code{\link{fitClere}}.} 
#' }
#' 
#' @seealso Overview : \code{\link{clere-package}} \cr 
#' Classes : \code{\linkS4class{Clere}} \cr 
#' Methods : \code{\link{plot}}, \code{\link{clusters}}, \code{\link{predict}}, \code{\link{summary}} \cr 
#' Functions : \code{\link{fitClere}} 
#' Datasets : \code{\link{numExpRealData}}, \code{\link{numExpSimData}}, \code{\link{algoComp}}
#' 
#' @export
#'
methods::setClass(
  Class = "Clere",
  representation = methods::representation(
    y = "numeric",
    x = "matrix",
    n = "integer",
    p = "integer",
    g = "numeric",
    nItMC = "numeric",
    nItEM = "numeric",
    nBurn = "numeric",
    dp = "numeric",
    nsamp = "numeric",
    sparse = "logical",
    analysis = "character",
    algorithm = "character",
    initialized = "logical",
    maxit = "numeric",
    tol = "numeric",
    seed = "integer",
    b = "numeric",
    pi = "numeric",
    sigma2 = "numeric",
    gamma2 = "numeric",
    intercept = "numeric",
    likelihood = "numeric",
    entropy = "numeric",
    P = "matrix",
    theta = "matrix",
    Zw = "matrix",
    Bw = "matrix",
    Z0 = "numeric",
    message = "character"
  ),
  prototype = methods::prototype(
    y = numeric(),
    x = matrix(NA, nrow = 0, ncol = 0),
    n = integer(),
    p = integer(),
    g = numeric(),
    nItMC = numeric(),
    nItEM = numeric(),
    nBurn = numeric(),
    dp = numeric(),
    nsamp = numeric(),
    sparse = logical(),
    analysis = character(),
    algorithm = character(),
    initialized = logical(),
    maxit = numeric(),
    tol = numeric(),
    seed = integer(),
    b = numeric(),
    pi = numeric(),
    sigma2 = numeric(),
    gamma2 = numeric(),
    intercept = numeric(),
    likelihood = numeric(),
    entropy = numeric(),
    P = matrix(NA, nrow = 0, ncol = 0),
    theta = matrix(NA, nrow = 0, ncol = 0),
    Zw = matrix(NA, nrow = 0, ncol = 0),
    Bw = matrix(NA, nrow = 0, ncol = 0),
    Z0 = numeric(),
    message = character()
  )
)


### Constructor ###
### Show ###
#' show
#' 
#' @name show
#' @aliases show,Clere-method
#' @docType methods
#' 
#' @rdname Clere-class
#'
#' @keywords internal
#' 
methods::setMethod(f = "show", signature = "Clere", definition = function(object) {
  showSlot <- function(slot) {
    sNames <- gsub("^[^@]*@(.*)", "\\1", slot)
    eSlot <- eval(parse(text = slot))
    tmp <- switch(EXPR = class(eSlot),
      "matrix" = {
        cat(paste0(" ~ ", sNames, " : [", nrow(eSlot), "x", ncol(eSlot), "]", collapse = ""))
        if (all(dim(eSlot) == 0)) {
          cat("NA")
        } else {
          cat("\n")
          nrowShow <- seq(min(5, nrow(eSlot)))
          ncolShow <- seq(min(5, ncol(eSlot)))
          shortObject <- eSlot[nrowShow, ncolShow]
          if (is.null(rownames(shortObject))) {
            rownames(shortObject) <- seq(nrow(shortObject))
          }
          if (is.null(colnames(shortObject))) {
            colnames(shortObject) <- seq(ncol(shortObject))
          }
          resFormat <- format(cbind(c("", rownames(shortObject)), rbind(colnames(shortObject), format(shortObject, digits = 4))), justify = "centre")
          if (nrow(shortObject) != nrow(eSlot)) {
            resFormat <- rbind(resFormat, c(".", sapply(seq(colnames(shortObject)), function(iCol) {
              paste0(rep(".", nchar(resFormat[1, 1])), collapse = "")
            })))
          }
          if (ncol(shortObject) != ncol(eSlot)) {
            resFormat <- cbind(resFormat, c(".", rep(paste0(rep(".", nchar(resFormat[1, 1])), collapse = ""), nrow(resFormat) - 1)))
          }
          cat(paste0("     ", apply(format(resFormat, justify = "centre"), 1, paste, collapse = " "), "\n", collapse = ""))
        }
        cat("\n")
      },
      "data.frame" = {
        cat(paste0(" ~ ", sNames, " : [", nrow(eSlot), "x", ncol(eSlot), "]", collapse = ""))
        if (all(dim(eSlot) == 0)) {
          cat("NA")
        } else {
          cat("\n")
          nrowShow <- seq(min(5, nrow(eSlot)))
          ncolShow <- seq(min(5, ncol(eSlot)))
          shortObject <- eSlot[nrowShow, ncolShow]
          if (is.null(rownames(shortObject))) {
            rownames(shortObject) <- seq(nrow(shortObject))
          }
          if (is.null(colnames(shortObject))) {
            colnames(shortObject) <- seq(ncol(shortObject))
          }
          resFormat <- format(cbind(c("", rownames(shortObject)), rbind(colnames(shortObject), format(shortObject, digits = 4))), justify = "centre")
          if (nrow(shortObject) != nrow(eSlot)) {
            resFormat <- rbind(resFormat, c(".", sapply(seq(colnames(shortObject)), function(iCol) {
              paste0(rep(".", nchar(resFormat[1, 1])), collapse = "")
            })))
          }
          if (ncol(shortObject) != ncol(eSlot)) {
            resFormat <- cbind(resFormat, c(".", rep(paste0(rep(".", nchar(resFormat[1, 1])), collapse = ""), nrow(resFormat) - 1)))
          }
          cat(paste0("     ", apply(format(resFormat, justify = "centre"), 1, paste, collapse = " "), "\n", collapse = ""))
        }
        cat("\n")
      },
      "numeric" = {
        cat(paste0(" ~ ", sNames, " : ", collapse = ""))
        if (length(eSlot) == 0) {
          cat("NA")
        } else {
          if (length(eSlot) > 1) {
            cat(paste0("[", length(eSlot), "] ", paste0(format(utils::head(eSlot), digits = 4), collapse = " ")))
          } else {
            cat(format(eSlot, digits = 4))
          }
        }
        cat("\n")
      },
      "character" = {
        cat(paste0(" ~ ", sNames, " : ", collapse = ""))
        if (length(eSlot) == 0) {
          cat("NA")
        } else {
          if (length(eSlot) > 1) {
            cat("[", length(eSlot), "] \"", paste0(utils::head(eSlot), collapse = "\" \""), "\"", sep = "")
          } else {
            cat(paste0("\"", eSlot, "\""))
          }
        }
        cat("\n")
      },
      {
        cat(paste0(" ~ ", sNames, " : ", collapse = ""))
        if (length(eSlot) == 0) {
          cat("NA")
        } else {
          if (length(eSlot) > 1) {
            cat(paste0("[", length(eSlot), "] ", paste0(utils::head(eSlot), collapse = " ")))
          } else {
            cat(eSlot)
          }
        }
        cat("\n")
      }
    )
    return(invisible())
  }
  showObject <- function(object) {
    cat("	~~~ Class:", class(object), "~~~\n")
    sNames <- paste0("object@", methods::slotNames(object))
    trash <- sapply(sNames, showSlot)
    return(invisible())
  }
  showObject(object)
  return(invisible(object))
})


#' Getteur
#'
#' @name [
#' @aliases [,Clere-method [,Clere,ANY,ANY,ANY-method
#' @docType methods
#' 
#' @rdname Clere-class
#'
#' @keywords internal
#' 
methods::setMethod(f = "[", signature = "Clere", definition = function(x, i, j, drop) {
  switch(EXPR = i,
    "y" = {
      if (missing(j)) {
        return(x@y)
      } else {
        if (j > length(x@y)) {
          stop("[Clere:get] indice out of limits", call. = FALSE)
        } else {
          return(x@y[j])
        }
      }
    },
    "x" = {
      return(x@x)
    },
    "n" = {
      if (missing(j)) {
        return(x@n)
      } else {
        if (j > length(x@n)) {
          stop("[Clere:get] indice out of limits", call. = FALSE)
        } else {
          return(x@n[j])
        }
      }
    },
    "p" = {
      if (missing(j)) {
        return(x@p)
      } else {
        if (j > length(x@p)) {
          stop("[Clere:get] indice out of limits", call. = FALSE)
        } else {
          return(x@p[j])
        }
      }
    },
    "g" = {
      if (missing(j)) {
        return(x@g)
      } else {
        if (j > length(x@g)) {
          stop("[Clere:get] indice out of limits", call. = FALSE)
        } else {
          return(x@g[j])
        }
      }
    },
    "nItMC" = {
      if (missing(j)) {
        return(x@nItMC)
      } else {
        if (j > length(x@nItMC)) {
          stop("[Clere:get] indice out of limits", call. = FALSE)
        } else {
          return(x@nItMC[j])
        }
      }
    },
    "nItEM" = {
      if (missing(j)) {
        return(x@nItEM)
      } else {
        if (j > length(x@nItEM)) {
          stop("[Clere:get] indice out of limits", call. = FALSE)
        } else {
          return(x@nItEM[j])
        }
      }
    },
    "nBurn" = {
      if (missing(j)) {
        return(x@nBurn)
      } else {
        if (j > length(x@nBurn)) {
          stop("[Clere:get] indice out of limits", call. = FALSE)
        } else {
          return(x@nBurn[j])
        }
      }
    },
    "dp" = {
      if (missing(j)) {
        return(x@dp)
      } else {
        if (j > length(x@dp)) {
          stop("[Clere:get] indice out of limits", call. = FALSE)
        } else {
          return(x@dp[j])
        }
      }
    },
    "nsamp" = {
      if (missing(j)) {
        return(x@nsamp)
      } else {
        if (j > length(x@nsamp)) {
          stop("[Clere:get] indice out of limits", call. = FALSE)
        } else {
          return(x@nsamp[j])
        }
      }
    },
    "sparse" = {
      if (missing(j)) {
        return(x@sparse)
      } else {
        if (j > length(x@sparse)) {
          stop("[Clere:get] indice out of limits", call. = FALSE)
        } else {
          return(x@sparse[j])
        }
      }
    },
    "analysis" = {
      if (missing(j)) {
        return(x@analysis)
      } else {
        if (j > length(x@analysis)) {
          stop("[Clere:get] indice out of limits", call. = FALSE)
        } else {
          return(x@analysis[j])
        }
      }
    },
    "algorithm" = {
      if (missing(j)) {
        return(x@algorithm)
      } else {
        if (j > length(x@algorithm)) {
          stop("[Clere:get] indice out of limits", call. = FALSE)
        } else {
          return(x@algorithm[j])
        }
      }
    },
    "initialized" = {
      if (missing(j)) {
        return(x@initialized)
      } else {
        if (j > length(x@initialized)) {
          stop("[Clere:get] indice out of limits", call. = FALSE)
        } else {
          return(x@initialized[j])
        }
      }
    },
    "maxit" = {
      if (missing(j)) {
        return(x@maxit)
      } else {
        if (j > length(x@maxit)) {
          stop("[Clere:get] indice out of limits", call. = FALSE)
        } else {
          return(x@maxit[j])
        }
      }
    },
    "tol" = {
      if (missing(j)) {
        return(x@tol)
      } else {
        if (j > length(x@tol)) {
          stop("[Clere:get] indice out of limits", call. = FALSE)
        } else {
          return(x@tol[j])
        }
      }
    },
    "seed" = {
      if (missing(j)) {
        return(x@seed)
      } else {
        if (j > length(x@seed)) {
          stop("[Clere:get] indice out of limits", call. = FALSE)
        } else {
          return(x@seed[j])
        }
      }
    },
    "b" = {
      if (missing(j)) {
        return(x@b)
      } else {
        if (j > length(x@b)) {
          stop("[Clere:get] indice out of limits", call. = FALSE)
        } else {
          return(x@b[j])
        }
      }
    },
    "pi" = {
      if (missing(j)) {
        return(x@pi)
      } else {
        if (j > length(x@pi)) {
          stop("[Clere:get] indice out of limits", call. = FALSE)
        } else {
          return(x@pi[j])
        }
      }
    },
    "sigma2" = {
      if (missing(j)) {
        return(x@sigma2)
      } else {
        if (j > length(x@sigma2)) {
          stop("[Clere:get] indice out of limits", call. = FALSE)
        } else {
          return(x@sigma2[j])
        }
      }
    },
    "gamma2" = {
      if (missing(j)) {
        return(x@gamma2)
      } else {
        if (j > length(x@gamma2)) {
          stop("[Clere:get] indice out of limits", call. = FALSE)
        } else {
          return(x@gamma2[j])
        }
      }
    },
    "intercept" = {
      if (missing(j)) {
        return(x@intercept)
      } else {
        if (j > length(x@intercept)) {
          stop("[Clere:get] indice out of limits", call. = FALSE)
        } else {
          return(x@intercept[j])
        }
      }
    },
    "likelihood" = {
      if (missing(j)) {
        return(x@likelihood)
      } else {
        if (j > length(x@likelihood)) {
          stop("[Clere:get] indice out of limits", call. = FALSE)
        } else {
          return(x@likelihood[j])
        }
      }
    },
    "entropy" = {
      if (missing(j)) {
        return(x@entropy)
      } else {
        if (j > length(x@entropy)) {
          stop("[Clere:get] indice out of limits", call. = FALSE)
        } else {
          return(x@entropy[j])
        }
      }
    },
    "P" = {
      return(x@P)
    },
    "theta" = {
      return(x@theta)
    },
    "Zw" = {
      return(x@Zw)
    },
    "Bw" = {
      return(x@Bw)
    },
    "Z0" = {
      if (missing(j)) {
        return(x@Z0)
      } else {
        if (j > length(x@Z0)) {
          stop("[Clere:get] indice out of limits", call. = FALSE)
        } else {
          return(x@Z0[j])
        }
      }
    },
    "message" = {
      if (missing(j)) {
        return(x@message)
      } else {
        if (j > length(x@message)) {
          stop("[Clere:get] indice out of limits", call. = FALSE)
        } else {
          return(x@message[j])
        }
      }
    },
    stop("[Clere:get] ", i, " is not a \"Clere\" slot", call. = FALSE)
  )
})


#' clusters method
#' 
#' This function makes returns the estimated clustering of variables.
#' 
#' 
#' @name clusters
#' @aliases clusters clusters-methods clusters,Clere-method
#' @docType methods
#' 
#' @param object [Clere]: Output object from \code{\link{fitClere}}.
#' @param threshold [numeric]: A numerical \code{threshold > 0.5}. If
#' \code{threshold = NULL} then the each variable is assigned to the cluster
#' having the largest associated posterior probability.
#' @param ... Additional arguments, not to be supplied in this version.
#' 
#' @seealso Overview : \code{\link{clere-package}} \cr 
#' Classes : \code{\linkS4class{Clere}} \cr 
#' Methods : \code{\link{show}}, \code{\link{plot}}, \code{\link{clusters}}, \code{\link{predict}}, \code{\link{summary}} \cr 
#' Functions : \code{\link{fitClere}} 
#' Datasets : \code{\link{numExpRealData}}, \code{\link{numExpSimData}}
#'
#' @export
#'
methods::setGeneric(name = "clusters", def = function(object, threshold = NULL, ...) {
  standardGeneric("clusters")
})
#' @export
methods::setMethod(f = "clusters", signature = "Clere", definition = function(object, threshold = NULL, ...) {
  if (is.null(threshold)) {
    return(apply(object@P, 1, which.max))
  } else {
    return(sapply(1:object@p, function(j) {
      counts <- sum(object@P[j, ] >= threshold)
      if (counts > 0) {
        if (counts > 1) {
          warning("[Clere:clusters] Variable ", j, " could have been assigned to other groups using ", threshold, " as threshold.\nYou may consider using a larger threshold or set option \" threshold = NULL\".", call. = TRUE)
        }
        return(which.max(object@P[j, ])[1])
      } else {
        return(NA)
      }
    }))
    ## return(apply(object@P, 1, function(x) ifelse(sum(x>=threshold)>0,which.max(x),NA)))
  }
})



#' predict method
#' 
#' This function makes prediction using a fitted model and a new matrix of
#' design. It returns a vector of predicted values of size equal to the number
#' of rows of matrix newx.
#' 
#' 
#' @name predict
#' @aliases predict predict-methods predict,Clere-method
#' @docType methods
#' 
#' @param object [Clere]: Output object from \code{\link{fitClere}}.
#' @param newx [matrix]: A numeric design matrix.
#' @param ... Additional arguments, not to be supplied in this version.
#' 
#' @seealso Overview : \code{\link{clere-package}} \cr 
#' Classes : \code{\linkS4class{Clere}} \cr 
#' Methods : \code{\link{show}}, \code{\link{plot}}, \code{\link{clusters}}, \code{\link{predict}}, \code{\link{summary}} \cr 
#' Functions : \code{\link{fitClere}} 
#' Datasets : \code{\link{numExpRealData}}, \code{\link{numExpSimData}}
#'
#' @export
#'
methods::setMethod(f = "predict", signature = "Clere", definition = function(object, newx, ...) {
  if (inherits(newx, "matrix")) {
    if (ncol(newx) == object@p) {
      return(object@intercept + rowMeans(newx %*% object@Bw, na.rm = TRUE))
    } else {
      stop(paste0("[Clere:predict] \"newx\" must be a matrix with ", object@p, " columns"), call. = FALSE)
    }
  } else {
    stop("[Clere:predict] \"newx\" is not a matrix", call. = FALSE)
  }
})


#' summary method
#' 
#' This function summarizes the output of function \code{\link{fitClere}}.
#' 
#' @name summary
#' @aliases summary summary-methods summary,Clere-method
#' @docType methods
#' 
#' @param object [Clere]: Output object from \code{\link{fitClere}}.
#' @param ... Additional arguments, not to be supplied in this version.
#' 
#' @seealso Overview : \code{\link{clere-package}} \cr 
#' Classes : \code{\linkS4class{Clere}} \cr 
#' Methods : \code{\link{show}}, \code{\link{plot}}, \code{\link{clusters}}, \code{\link{predict}}, \code{\link{summary}} \cr 
#' Functions : \code{\link{fitClere}} 
#' Datasets : \code{\link{numExpRealData}}, \code{\link{numExpSimData}}
#'
#' @export
#'
methods::setMethod(f = "summary", signature = "Clere", definition = function(object, ...) {
  if (missing(object)) {
    stop("[Clere:summary] \"object\" is missing", call. = FALSE)
  }
  nbVarGroups <- length(object@b)
  K <- 2 * (nbVarGroups + 1)
  nd <- 4
  sep <- "\t"
  summaryClere <- methods::new("sClere",
    analysis = object@analysis,
    g = object@g,
    nbVarGroups = length(object@b),
    algorithm = object@algorithm,
    intercept = object@intercept,
    b = object@b,
    pi = object@pi,
    sigma2 = object@sigma2,
    gamma2 = object@gamma2,
    likelihood = object@likelihood,
    entropy = object@entropy,
    AIC = -2 * object@likelihood + 2 * K,
    BIC = -2 * object@likelihood + K * log(object@n),
    ICL = -2 * object@likelihood + K * log(object@n) + object@entropy
  )
  return(summaryClere)
})


#' plot method
#' 
#' Graphical summary for MCEM/SEM-Gibbs estimation.  This function represents
#' the course of the model parameters in view of the iterations of the
#' estimation algorithms implemented in \code{\link{fitClere}}.
#' 
#' 
#' @name plot-methods
#' @aliases plot plot-methods plot,Clere-method plot,Clere,ANY-method
#' @docType methods
#' 
#' @param x [Clere]: Output object from \code{\link{fitClere}}.
#' @param y [any]: Unused parameter.
#' @param ... Additional arguments, not to be supplied in this version.
#' 
#' @seealso Overview : \code{\link{clere-package}} \cr 
#' Classes : \code{\linkS4class{Clere}}, \code{\linkS4class{Pacs}} \cr 
#' Methods : \code{\link{plot}}, \code{\link{clusters}}, \code{\link{predict}}, \code{\link{summary}} \cr 
#' Functions : \code{\link{fitClere}}, \code{\link{fitPacs}} 
#' Datasets : \code{\link{numExpRealData}}, \code{\link{numExpSimData}}, \code{\link{algoComp}}
#' 
#' @export
#' 
# setGeneric(name = "plot", def = function(x, ...) {standardGeneric("plot")})
methods::setMethod(f = "plot", signature = "Clere", definition = function(x, ...) {
  if (nrow(x@theta) >= 2) {
    snbRow <- seq(nrow(x@theta))
    K <- length(x@b)
    sK <- seq(K)
    op <- graphics::par(mfrow = c(2, 2))
    xlabname <- paste(x@algorithm, "iterations")

    ### First plot
    graphics::matplot(snbRow, x@theta[, 1 + (sK)], type = "l", ylab = "The b's", lty = 1, xlab = xlabname)
    graphics::abline(v = x@nBurn, col = "grey")

    ### Second plot
    graphics::matplot(snbRow, x@theta[, 1 + K + (sK)], type = "l", ylab = "The pi's", lty = 1, xlab = xlabname)
    graphics::abline(v = x@nBurn, col = "grey")

    ### Third plot
    graphics::matplot(snbRow, x@theta[, c(2 * K + 2, 2 * K + 3)],
      type = "l", ylab = "sigma^2 and gamma^2", lty = 1, xlab = xlabname,
      ylim = c(0.5 * min(x@theta[, c(2 * K + 2, 2 * K + 3)], na.rm = TRUE), 1.5 * max(x@theta[, c(2 * K + 2, 2 * K + 3)], na.rm = TRUE))
    )
    graphics::legend("topright", legend = c("sigma2", "gamma2"), box.lty = 0, lty = 1, col = 1:2)
    graphics::abline(v = x@nBurn, col = "grey")

    ### Fourth plot
    graphics::plot(snbRow, x@theta[, 1], type = "l", ylab = "Intercept", lty = 1, xlab = xlabname)
    graphics::abline(v = x@nBurn, col = "grey")

    graphics::par(op)
  } else {
    stop("[Clere:plot] SEM iterations is below 2", call. = FALSE)
  }
})
