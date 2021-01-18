### var2dummy.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: jun 22 2017 (16:03) 
## Version: 
## last-updated: sep 28 2018 (11:57) 
##           By: Brice Ozenne
##     Update #: 41
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * var2dummy
#' @title Convert Variable Names to Dummy Variables Names.
#' @description When dealing with categorical variables, the \code{estimate} function convert the categorical variables into dummy variables.
#' This function convert a set of variable names to their corresponding name in the model with dummy variables
#' @name var2dummy
#' 
#' @param object a \code{lvm} object.
#' @param var [character] the variable to be transformed.
#' @param data [data.frame] dataset according to which the model should be updated.
#' @param rm.first.factor [logical] should the first level of each categorical variable be ignored?
#' @param ... [internal] additional arguments to be passed from \code{var2dummy.lvm} to \code{var2dummy.list}.
#' 
#' @examples
#'
#' \dontrun{
#' var2dummy <- lavaSearch2:::var2dummy
#' var2dummy.list <- lavaSearch2:::var2dummy.list
#' var2dummy.lvm <- lavaSearch2:::var2dummy.lvm
#' 
#' m <- lvm()
#' regression(m) <- c(y1,y2,y3)~u
#' regression(m) <- u ~ X1+X2
#' var2dummy(m, var = c("X1","X2"))
#' categorical(m,labels=c("M","F","MF")) <- ~X1
#' var2dummy(m, var = c("X1","X2"))
#' categorical(m,labels=c("1","2","3")) <- ~X2
#' var2dummy(m, var = c("X1","X2"))
#' }
#' 
#' @keywords internal
`var2dummy` <-
    function(object,...) UseMethod("var2dummy")

## * var2dummy.list
#' @rdname var2dummy
#' @export
var2dummy.list <- function(object, var, rm.first.factor = TRUE, ...){

    var <- stats::setNames(var,var)
    ## convertion to dummy variable name for categorical variables
    factor.var <- names(object$x$attributes$labels)
    
    if(!is.null(var) && any(var %in% factor.var)){
        subvar <- var[var %in% factor.var]
        for(iFactor in subvar){ # iFactor <- "X1"
            newvar <- paste0(iFactor,object$x$attributes$labels[[iFactor]])
            if(rm.first.factor){newvar <- newvar[-1]}
            newvar <- stats::setNames(newvar, rep(iFactor, length(newvar)))
            var <- c(var[names(var)!=iFactor],newvar)            
        }
    }
    return(var)
}

## * var2dummy.lvm
#' @rdname var2dummy
#' @export
var2dummy.lvm <- function(object, data = NULL, ...){

    if(is.null(data)){
        data <- lava::sim(object, n = 1)
    }
    object2 <- lava_categorical2dummy(object, data)
    
    ## recover attributes for models not defined using categorical
    obsvars <- setdiff(lava::vars(object2$x),lava::latent(object))
    
    if(any(obsvars %in% names(data) == FALSE)){
        missing.vars <- obsvars[obsvars %in% names(data) == FALSE]
        test.num <- sapply(1:NCOL(data), function(col){is.numeric(data[[col]])})
        possible.match <- names(data)[test.num==FALSE]
        n.possible.match <- length(possible.match)
        
        ls.labels <- list()
        for(iMatch in 1:n.possible.match){ # iMatch <- 1
            iVar <- possible.match[iMatch]
            iLabel <- levels(as.factor(data[[iVar]]))
            if(any(paste0(iVar,iLabel) %in% missing.vars)){
                ls.labels[[iVar]] <- iLabel
            }
        }

        object2$x$attributes$labels <- ls.labels
    }
    
    res <- var2dummy(object2, ...)
    return(res)
}

#----------------------------------------------------------------------
### var2dummy.R ends here
