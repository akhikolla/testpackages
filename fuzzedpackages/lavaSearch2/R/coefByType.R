### coefByType.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb  5 2018 (11:37) 
## Version: 
## Last-Updated: okt  4 2018 (17:42) 
##           By: Brice Ozenne
##     Update #: 32
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation
#' @title Extract the Coefficient by Type
#' @description Extract specific types of coefficient from a \code{lvm} object:
#' covariance coefficient(s) (\code{coefCov}),
#' extra parameter(s) (\code{coefExtra}),
#' position in the list of models for each coefficient (\code{coefIndexModel}),
#' intercept coefficient(s) (\code{coefIntercept}),
#' coefficient(s) that are used as reference (\code{coefRef}),
#' regression coefficient(s) (\code{coefReg}),
#' variance coefficient(s) (\code{coefVar}).
#'  
#' @name coefByType
#' 
#' @param object a lvm model or a fitted lvm model 
#' @param value should the name of the coefficient be returned? Else return the coefficients
#' @param keep.var should the variance coefficients be returned?
#' @param ... arguments to be passed to \code{}
#'
#' @return A vector containing the names of the positions of the coefficients.
#' 
#' @examples 
#' #### regression ####
#' m <- lvm(Y~X1+X2)
#' e <- estimate(m, lava::sim(m, 1e2))
#' 
#' coefCov(m)
#' coefCov(m, value = TRUE)
#'
#' coefCov(m, keep.var = TRUE)
#' coefCov(m, value = TRUE, keep.var = TRUE)
#'
#' coefIndexModel(m)
#' coefIndexModel(e)
#' 
#' coefIntercept(m)
#' coefIntercept(m, value = TRUE)
#'
#' coefReg(m)
#' coefReg(m, value = TRUE)
#' 
#' #### LVM ####
#' m <- lvm()
#' regression(m) <- c(y1,y2,y3)~u
#' regression(m) <- u~x1+x2
#' latent(m) <- ~u
#' covariance(m) <- y1~y2
#'
#' m.Sim <- m
#' categorical(m.Sim, labels = c("a","b","c")) <- ~x2
#' e <- estimate(m, lava::sim(m.Sim, 1e2))
#'
#' coefCov(m)
#' coefCov(m, value = TRUE) 
#'
#' coefCov(m, keep.var = TRUE)
#' coefCov(m, value = TRUE, keep.var = TRUE)
#' 
#' coefExtra(m)
#'
#' coefIndexModel(m)
#' coefIndexModel(e)
#'
#' ## additional categorical variable 
#' categorical(m, labels = as.character(1:3)) <- "X1"
#' 
#' coefExtra(m)
#' coefExtra(m, value = TRUE)
#'
#' ## additional categorical variable
#' categorical(m, labels = as.character(1:3)) <- "x1"
#'
#' coefIntercept(m)
#' coefIntercept(m, value = TRUE)
#' coefIntercept(e)
#'
#' coefReg(e, value = TRUE)
#' 
#' #### multigroup ####
#' m <- lvm(Y~X1+X2)
#' eG <- estimate(list(m,m), list(lava::sim(m, 1e2), lava::sim(m, 1e2)))
#' 
#' coefIndexModel(eG)
#'
#' @concept extractor

## * coefCov
#' @rdname coefByType
#' @export
`coefCov` <-
  function(object, value, keep.var, ...) UseMethod("coefCov")

## ** coefCov.lvm
#' @rdname coefByType
#' @export
coefCov.lvm <- function(object, value = FALSE, keep.var = FALSE, ...){

    res <- retainType(type = coefType(object, ...),
                      validType = c("covariance", if(keep.var){"variance"}else{NULL}),
                      value = value)

    return(res)
}

## ** coefCov.lvmfit
#' @rdname coefByType
#' @export
coefCov.lvmfit <- coefCov.lvm

## ** coefCov.multigroup
#' @rdname coefByType
#' @export
coefCov.multigroup <- coefCov.lvm

## * coefExtra

#' @rdname coefByType
#' @export
`coefExtra` <-
  function(object, value, ...) UseMethod("coefExtra")

## ** coefExtra.lvm
#' @rdname coefByType
#' @export
coefExtra.lvm <- function(object, value = FALSE, ...){
    
    res <- retainType(type = coefType(object, ...),
                      validType = "extra",
                      value = value) 
    
    return(res)    
}

## ** coefExtra.lvmfit
#' @rdname coefByType
#' @export
coefExtra.lvmfit <- coefExtra.lvm

## ** coefExtra.multigroup
#' @rdname coefByType
#' @export
coefExtra.multigroup <- coefExtra.lvm

## * coefIndexModel
#' @rdname coefByType
#' @export
`coefIndexModel` <-
  function(object,...) UseMethod("coefIndexModel")

## ** coefIndexModel.lvm
#' @rdname coefByType
#' @export
coefIndexModel.lvm <- function(object, ...){
    name.coef <- stats::coef(object)
    index <- rep(1, length(name.coef))
    names(index) <- name.coef
    return(index)
}

## ** coefIndexModel.lvmfit
#' @rdname coefByType
#' @export
coefIndexModel.lvmfit <- function(object, ...){
    name.coef <- names(stats::coef(object))
    index <- rep(1, length(name.coef))
    names(index) <- name.coef
    return(index)
}

## ** coefIndexModel.multigroup
#' @rdname coefByType
#' @export
coefIndexModel.multigroup <- function(object, ...){
    n.model <- length(object$lvm)

    ## new coef names
    allCoef <- object$name
    n.allCoef <- length(allCoef)
    index.AllCoef <- object$coef
  
    index <- stats::setNames(rep(NA, n.allCoef), allCoef)
  
    for(iModel in 1:n.model){ # iModel <- 1
        index[index.AllCoef[[iModel]]] <- iModel
    }
  
    #### export
    return(index)
}

## ** coefIndexModel.multigroupfit
#' @rdname coefByType
#' @export
coefIndexModel.multigroupfit <- function(object, ...){
    out <- coefIndexModel(object$model)
    return(out[names(stats::coef(object))])
}
  
## * coefIntercept

#' @rdname coefByType
#' @export
`coefIntercept` <-
  function(object, value, ...) UseMethod("coefIntercept")

## ** coefIntercept.lvm
#' @rdname coefByType
#' @export
coefIntercept.lvm <- function(object, value = FALSE, ...){ 

    res <- retainType(type = coefType(object, ...),
                      validType = "intercept",
                      value = value)

    return(res)
}

## ** coefIntercept.lvmfit
#' @rdname coefByType
#' @export
coefIntercept.lvmfit <- coefIntercept.lvm

## ** coefIntercept.multigroup
#' @rdname coefByType
#' @export
coefIntercept.multigroup <- coefIntercept.lvm

## * coefRef
#' @rdname coefByType
#' @export
`coefRef` <-
  function(object, value, ...) UseMethod("coefRef")

## ** coefRef.lvmfit
#' @rdname coefByType
#' @export
coefRef.lvmfit <- function(object, value = FALSE, ...){
    
    res <- retainType(type = attr(coefType(object, ...), "reference"),
                      validType = TRUE,
                      value = value)

    return(res)    
}

## * coefReg
#' @rdname coefByType
#' @export
`coefReg` <-
  function(object, value, ...) UseMethod("coefReg")

## ** coefReg.lvm
#' @rdname coefByType
#' @export
coefReg.lvm <- function(object, value = FALSE, ...){
    
     res <- retainType(type = coefType(object, ...),
                      validType = "regression",
                      value = value)

     return(res)
}

## ** coefReg.lvmfit
#' @rdname coefByType
#' @export
coefReg.lvmfit <- coefReg.lvm

## ** coefReg.multigroup
#' @rdname coefByType
#' @export
coefReg.multigroup <- coefReg.lvm

## * coefVar

#' @rdname coefByType
#' @export
`coefVar` <-
  function(object, value, ...) UseMethod("coefVar")

## ** coefVar.lvm
#' @rdname coefByType
#' @export
coefVar.lvm <- function(object, value = FALSE, ...){ 

    res <- retainType(type = coefType(object, ...),
                      validType = "variance",
                      value = value)

    return(res)
}

## ** coefVar.lvmfit
#' @rdname coefByType
#' @export
coefVar.lvmfit <- coefVar.lvm

## ** coefVar.multigroup
#' @rdname coefByType
#' @export
coefVar.multigroup <- coefVar.lvm

## * retainType  (needed for coefCov/Latent/Ref)
retainType <- function(type, validType, value){
  index.var <- which(type %in% validType)
  
  if(length(index.var)>0){
      if(value){
          return(names(type)[index.var])
      }else{
          return(index.var)
      }
  }else{
      return(NULL)
  }
}


##----------------------------------------------------------------------
### coefByType.R ends here
