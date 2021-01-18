### conditionalMoment.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 27 2017 (16:59) 
## Version: 
## last-updated: feb  8 2019 (11:47) 
##           By: Brice Ozenne
##     Update #: 1139
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * conditionalMoment - documentation
#' @title Prepare the Computation of score2
#' @description Compute the conditional mean and variance,
#' and their first and second derivative regarding the model parameters.
#' @name conditionalMoment
#' 
#' @param object,x a latent variable model.
#' @param data [data.frame] data set.
#' @param formula [formula] two-sided linear formula.
#' @param param,p [numeric vector] the fitted coefficients.
#' @param attr.param [character vector] the type of each coefficient
#' (e.g. mean or variance coefficient).
#' @param ref.group [character vector] the levels of the variable defining the variance component in a generic covariance matrix.
#' @param second.order [logical] should the terms relative to the third derivative of the likelihood be be pre-computed?
#' @param cluster [integer vector] the grouping variable relative to which the observations are iid.
#' @param n.cluster [integer >0] the number of i.i.d. observations.
#' @param n.endogenous [integer >0] the number of outcomes.
#' @param usefit,value [logical] If TRUE the coefficients estimated by the model are used to pre-compute quantities. Only for lvmfit objects.
#' @param name.endogenous [character vector, optional] name of the endogenous variables
#' @param name.latent [character vector, optional] name of the latent variables
#' @param ... [internal] only used by the generic method or by the <- methods.
#' 
#' @details For lvmfit objects, there are two levels of pre-computation:
#' \itemize{
#' \item a basic one that do no involve the model coefficient (\code{conditionalMoment.lvm}).
#' \item an advanced one that require the model coefficients (\code{conditionalMoment.lvmfit}). 
#' }
#' 
#' @examples
#' m <- lvm(Y1~eta,Y2~eta,Y3~eta)
#' latent(m) <- ~eta
#'
#' d <- lava::sim(m,1e2)
#' e <- estimate(m, d)
#'
#' ## basic pre-computation
#' res1 <- conditionalMoment(e, data = d,
#'                          first.order = FALSE, second.order = FALSE,
#'                          name.endogenous = endogenous(e),
#'                          name.latent = latent(e), usefit = FALSE)
#' res1$skeleton$Sigma
#' 
#' ## full pre-computation
#' res2 <- conditionalMoment(e, param = coef(e), data = d,
#'                          first.order = FALSE, second.order = FALSE,
#'                          name.endogenous = endogenous(e),
#'                          name.latent = latent(e), usefit = TRUE
#' )
#' res2$value$Sigma
#'
#' @concept small sample inference
#' @concept derivative of the score equation
#' 
#' @keywords internal
#' @export
`conditionalMoment` <-
  function(object, ...) UseMethod("conditionalMoment")


## * conditionalMoment.lm
#' @rdname conditionalMoment
#' @export
conditionalMoment.lm <- function(object, data, param,
                                 name.endogenous,
                                 first.order, second.order, ...){

    out <- list(param = param,
                name.3deriv = "sigma2")
    
    ## design matrix
    X <- model.matrix(formula(object), data)

    ## linear predictor
    out$mu <- X %*% param[colnames(X)]
    
    ## residuals variance
    out$Omega <- matrix(param["sigma2"], nrow = 1, ncol = 1,
                    dimnames = list(name.endogenous, name.endogenous))
    
    ## ** first order
    if(first.order){
        out$dmu <- lapply(1:NCOL(X), function(i){
            M <- X[,i,drop=FALSE]
            colnames(M) <- name.endogenous
            return(M)
        })
        names(out$dmu) <- colnames(X)
        out$dOmega = list(sigma2 = matrix(1))
    }

    ## ** second order
    if(second.order){
        out$d2mu <- NULL
        out$d2Omega <- NULL
    }

    ## ** export
    return(out)
}

## * conditionalMoment.gls
#' @rdname conditionalMoment
#' @export
conditionalMoment.gls <- function(object, data, formula, 
                                  param, attr.param, ref.group,
                                  first.order, second.order,
                                  index.Omega, vec.OmegaMat, cluster, n.cluster,
                                  name.endogenous, n.endogenous,
                                  ...){

    if(first.order == FALSE && second.order == TRUE){
        stop("Cannot pre-compute quantities for the second order derivatives ",
             "without those for the first order derivatives \n")
    }
### ** prepare

    ## *** coefficients
    name.varcoef <- attr.param$var.coef
    name.corcoef <- attr.param$cor.coef
    n.varcoef <- length(name.varcoef)
    n.corcoef <- length(name.corcoef)
    var.coef <- param[name.varcoef]
    cor.coef <- param[name.corcoef]
        
    class.var <- class(object$modelStruct$varStruct)
    class.cor <- class(object$modelStruct$corStruct)

    
    ## *** design matrix    
    X <- stats::model.matrix(formula, data)
    X <- X[,attr.param$mean.coef,drop=FALSE] ## drop unused columns (e.g. factor with 0 occurence)    
    attr(X,"assign") <- NULL
    attr(X,"contrasts") <- NULL
    
    ## *** variance terms
    if("NULL" %in% class.var == FALSE){
        name.otherVar <- setdiff(names(var.coef),"sigma2")
        factor.varcoef <- setNames(c(1,var.coef[name.otherVar]),
                                   attr(object$modelStruct$varStruct,"groupNames"))
        sigma2.base0 <- factor.varcoef[ref.group]        
    }else{
        name.otherVar <- NULL
        sigma2.base0 <- setNames(rep(1, n.endogenous), name.endogenous)
    }
    sigma2.base <- sigma2.base0 * var.coef["sigma2"]

    ## *** corelation terms
    if("NULL" %in% class.cor == FALSE){
        M.corcoef <- matrix("", n.endogenous, n.endogenous,
                            dimnames = list(name.endogenous,name.endogenous))
        M.corcoef[which(lower.tri(M.corcoef))] <- name.corcoef
        M.corcoef <- symmetrize(M.corcoef)

        index.lower.tri <- which(lower.tri(M.corcoef))
        indexArr.lower.tri <- which(lower.tri(M.corcoef), arr.ind = TRUE)

        Msigma2.base0 <- matrix(0, n.endogenous, n.endogenous,
                                dimnames = list(name.endogenous, name.endogenous))
        Msigma2.base0[index.lower.tri] <- apply(indexArr.lower.tri, 1, function(x){sqrt(prod(sigma2.base0[x]))})
        Msigma2.base0 <- symmetrize(Msigma2.base0)
    }else{
        M.corcoef <- NULL
        Msigma2.base0 <- NULL
        index.lower.tri <- NULL
        indexArr.lower.tri <- NULL
    }

    ## *** export
    out <- list(param = param,
                name.3deriv = c(name.varcoef,name.corcoef),
                skeleton = list(class.cor = class.cor,
                                class.var = class.var,
                                sigma2.base0 = sigma2.base0,
                                Msigma2.base0 = Msigma2.base0,
                                ref.group = ref.group,
                                n.endogenous = n.endogenous,
                                name.endogenous = name.endogenous,
                                M.corcoef = M.corcoef,
                                index.lower.tri = index.lower.tri,
                                indexArr.lower.tri = indexArr.lower.tri,
                                var.coef = var.coef,
                                name.varcoef = name.varcoef,
                                n.varcoef = n.varcoef,
                                name.otherVar = name.otherVar,
                                ref.group = ref.group,
                                cor.coef = cor.coef,
                                name.corcoef = name.corcoef,
                                n.corcoef = n.corcoef,
                                cluster  = cluster,
                                n.cluster = n.cluster))
    
### ** Reconstruct conditional mean
    ## transpose necessary because of the way index.OmegaMat was computed
    out$mu <- matrix(NA, nrow = n.cluster, ncol = n.endogenous,
                     dimnames = list(NULL, name.endogenous))
    out$mu[vec.OmegaMat] <- X %*% param[colnames(X)]
    
    
### ** Reconstruct conditional variance covariance matrix
    out$Omega <- .getVarCov2(object,
                             param = param,
                             attr.param = attr.param,
                             name.endogenous = name.endogenous,
                             n.endogenous = n.endogenous,
                             ref.group = ref.group)
    
### ** first order
    if(first.order){
        outD1 <- skeletonDtheta(object,
                                class.cor = class.cor,
                                class.var = class.var,
                                X = X,  
                                sigma2.base0 = sigma2.base0,
                                Msigma2.base0 = Msigma2.base0,
                                M.corcoef = M.corcoef,
                                ref.group = ref.group,
                                name.endogenous = name.endogenous,
                                n.endogenous = n.endogenous,
                                index.lower.tri = index.lower.tri,
                                indexArr.lower.tri = indexArr.lower.tri,
                                cluster = cluster,
                                n.cluster = n.cluster,
                                var.coef = var.coef,
                                name.varcoef = name.varcoef,
                                name.otherVar = name.otherVar,
                                n.varcoef = n.varcoef,
                                cor.coef = cor.coef,
                                name.corcoef = name.corcoef,
                                n.corcoef = n.corcoef,
                                index.Omega = index.Omega,
                                update.mean = TRUE, update.variance = TRUE,
                                ...) ## ... to pass coef.rancoef
        out$dmu <- outD1$dmu
        out$dOmega <- outD1$dOmega
    }
    
### ** second order
    if(second.order){
        out$d2mu <- NULL
        out$d2Omega <- skeletonDtheta2(object,
                                       dOmega = out$dOmega,
                                       class.cor = class.cor,
                                       class.var = class.var,
                                       M.corcoef = M.corcoef,
                                       n.endogenous = n.endogenous,
                                       index.lower.tri = index.lower.tri,
                                       indexArr.lower.tri = indexArr.lower.tri,
                                       var.coef = var.coef,
                                       name.varcoef = name.varcoef,
                                       name.otherVar = name.otherVar,
                                       n.varcoef = n.varcoef,
                                       cor.coef = cor.coef,
                                       name.corcoef = name.corcoef)
    }
    
### ** export
    return(out)
    
}

## * conditionalMoment.lme
#' @rdname conditionalMoment
#' @export
conditionalMoment.lme <- function(object, attr.param, ...){

    name.rancoef <- attr.param$ran.coef
    out <- conditionalMoment.gls(object, attr.param = attr.param,
                                 name.rancoef = name.rancoef, ...)

    ##  the derivative regarding the random effect is added by skeletonDtheta.lme
    out$name.3deriv <- c(out$name.3deriv, name.rancoef)
    out$skeleton$class.ran <- class(object$modelStruct$reStruct)

    return(out)
}
## * conditionalMoment.lvm
#' @rdname conditionalMoment
#' @export
conditionalMoment.lvm <- function(object, data,
                                  first.order, second.order,
                                  name.endogenous, name.latent,
                                  ...){

    if(first.order == FALSE && second.order == TRUE){
        stop("Cannot pre-compute quantities for the second order derivatives ",
             "without those for the first order derivatives \n")
    }

### ** Initialize conditional moments   
    Moment <- skeleton(object,
                       name.endogenous = name.endogenous, 
                       name.latent = name.latent, 
                       as.lava = TRUE)

### ** Initialize partial derivatives of the conditional moments
    if(first.order){
        dMoment <- skeletonDtheta(object, data = data,
                                  df.param.all = Moment$df.param,
                                  param2originalLink = Moment$param2originalLink,
                                  name.endogenous = name.endogenous, 
                                  name.latent = name.latent)
    }else{
        dMoment <- NULL
    }
    
### ** Initialize second order partial derivatives of the conditional moments
    if(second.order){
        d2Moment <- skeletonDtheta2(object, data = data,
                                    df.param.all = Moment$df.param,
                                    param2originalLink = Moment$param2originalLink,
                                    name.latent = name.latent)
    }else{
        d2Moment <- NULL
    }

### ** Export
    return(c(Moment, list(dMoment.init = dMoment, d2Moment.init = d2Moment)))
}
    
    
## * conditionalMoment.lvmfit
#' @rdname conditionalMoment
#' @export
conditionalMoment.lvmfit <- function(object, data, param, 
                                     first.order, second.order, usefit,
                                     ...){

### ** normalize arguments
    name.endogenous <- endogenous(object)
    n.endogenous <- length(name.endogenous)
    name.latent <- latent(object)
    n.latent <- length(name.latent)

    data <- as.matrix(data[,lava::manifest(object),drop=FALSE])

### ** initialize
    if(is.null(object$conditionalMoment)){       
        object$conditionalMoment <- conditionalMoment(lava::Model(object),
                                                      data = data,
                                                      first.order = first.order,
                                                      second.order = second.order,
                                                      name.endogenous = name.endogenous,
                                                      name.latent = name.latent)

        ##  param with non-zero third derivative
        type.3deriv <- c("alpha","Gamma","Lambda","B","Psi_var","Sigma_var","Psi_cov","Sigma_cov")
        index.keep <- intersect(which(!is.na(object$conditionalMoment$df.param$lava)),
                                which(object$conditionalMoment$df.param$detail %in% type.3deriv)
                                )    
        object$conditionalMoment$name.3deriv <- object$conditionalMoment$df.param[index.keep, "originalLink"]
    }

### ** update according to the value of the model coefficients
    if(usefit){

        ## *** conditional moments
        object$conditionalMoment$value <- skeleton(object, data = data, param = param,
                                                   name.endogenous = name.endogenous, 
                                                   name.latent = name.latent)

        if(object$conditionalMoment$skeleton$toUpdate["param"]){
            object$conditionalMoment$param <- coef(object)
        }
        if(object$conditionalMoment$skeleton$toUpdate["mu"]){            
            if(n.latent==0){
                object$conditionalMoment$mu <- object$conditionalMoment$value$nu.XK
            }else{
                object$conditionalMoment$mu <- object$conditionalMoment$value$nu.XK + object$conditionalMoment$value$alpha.XGamma.iIB %*% object$conditionalMoment$value$Lambda
            }            
        }
        if(object$conditionalMoment$skeleton$toUpdate["Omega"]){
            object$conditionalMoment$Omega <- getVarCov2(object)
        }
        
        ## *** first order derivatives
        if(first.order){            
            out <- skeletonDtheta(object,
                                  name.endogenous = name.endogenous, 
                                  name.latent = name.latent)
            object$conditionalMoment$dmu <- out$dmu
            object$conditionalMoment$dOmega <- out$dOmega            
        }

        ## *** second order derivatives
        if(second.order){
            out2 <- skeletonDtheta2(object)
            object$conditionalMoment$d2mu <- out2$d2mu
            object$conditionalMoment$d2Omega <- out2$d2Omega
        }
       
    }
     
### ** Export
    return(object$conditionalMoment)
}

