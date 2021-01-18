### getVarCov2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 27 2018 (09:55) 
## Version: 
## Last-Updated: jul 25 2019 (10:13) 
##           By: Brice Ozenne
##     Update #: 33
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * getVarCov2

#' @title Reconstruct the Conditional Variance Covariance Matrix
#' @description Reconstruct the conditional variance covariance matrix from a nlme or lvm model.
#' Only compatible with specific correlation and variance structure.
#' @name getVarCov2
#'
#' @param object a \code{gls} or \code{lme} object
#' @param param [numeric vector] values for the model parameters.
#' @param data [data.frame] the data set.
#' @param cluster [integer vector] the grouping variable relative to which the observations are iid.
#' @param ... [internal] only used by the generic method.
#' 
#' @details The compound symmetry variance-covariance matrix in a gls model is of the form:
#' \tabular{cccc}{
#' \eqn{\Sigma =} \tab \eqn{\sigma^2} \tab \eqn{\sigma^2 \rho} \tab \eqn{\sigma^2 \rho} \cr
#' \tab . \tab \eqn{\sigma^2} \tab \eqn{\sigma^2 \rho} \cr
#' \tab . \tab . \tab \eqn{\sigma^2}
#' }
#'
#' The unstructured variance-covariance matrix in a gls model is of the form:
#'  \tabular{cccc}{
#' \eqn{\Sigma =} \tab \eqn{\sigma^2} \tab \eqn{\sigma^2 \sigma_2 \rho_{1,2}} \tab \eqn{\sigma^2 \sigma_3 \rho_{1,3}} \cr
#' \tab . \tab \eqn{\sigma^2 \sigma_2^2} \tab \eqn{\sigma^2 \sigma_2 \sigma_3 \rho_{2,3}} \cr
#' \tab . \tab . \tab \eqn{\sigma^2 \sigma_3^2}
#' }
#' @return A list containing the residual variance-covariance matrix in the element Omega.
#' 
#' @examples
#' 
#' ## simulate data 
#' library(nlme)
#' n <- 5e1
#' mSim <- lvm(c(Y1~1*eta,Y2~1*eta,Y3~1*eta,eta~G))
#' latent(mSim) <- ~eta
#' transform(mSim,Id~Y1) <- function(x){1:NROW(x)}
#' set.seed(10)
#' dW <- lava::sim(mSim,n,latent = FALSE)
#' dW <- dW[order(dW$Id),,drop=FALSE]
#' dL <- reshape2::melt(dW,id.vars = c("G","Id"), variable.name = "time")
#' dL <- dL[order(dL$Id),,drop=FALSE]
#' dL$Z1 <- rnorm(NROW(dL))
#' dL$time.num <- as.numeric(as.factor(dL$time))
#' 
#' #### iid model #### 
#' e1.gls <- nlme::gls(Y1 ~ G, data = dW, method = "ML")
#' getVarCov2(e1.gls, cluster = 1:n)$Omega
#' 
#' #### heteroschedasticity ####
#' dW$group <- rbinom(n, size = 1, prob = 1/2)
#' dW$repetition <- as.numeric(as.factor(dW$group))
#' e2a.gls <- nlme::gls(Y1 ~ G, data = dW, method = "ML",
#'                     weights = varIdent(form =~ repetition|group))
#' getVarCov2(e2a.gls, cluster = 1:n)$Omega
#'
#' 
#' e2b.gls <- nlme::gls(value ~ 0+time + time:G,
#'                    weight = varIdent(form = ~ time.num|time),
#'                    data = dL, method = "ML")
#' getVarCov2(e2b.gls, cluster = "Id")$Omega
#'
#' #### compound symmetry ####
#' e3.gls <- nlme::gls(value ~ time + G,
#'                    correlation = corCompSymm(form = ~1| Id),
#'                    data = dL, method = "ML")
#' getVarCov2(e3.gls)$Omega
#' 
#' #### unstructured ####
#' e4.gls <- nlme::gls(value ~ time,
#'                     correlation = corSymm(form = ~time.num| Id),
#'                     weight = varIdent(form = ~ 1|time),
#'                     data = dL, method = "ML")
#' getVarCov2(e4.gls)$Omega
#'
#' #### lvm model ####
#' m <- lvm(c(Y1~1*eta,Y2~1*eta,Y3~1*eta,eta~G))
#' latent(m) <- ~eta
#' e <- estimate(m, dW)
#' getVarCov2(e)
#' 
#' @concept extractor
#' 
#' @export
`getVarCov2` <-
    function(object, ...) UseMethod("getVarCov2")

## * getVarCov2.gls
#' @rdname getVarCov2
#' @export
getVarCov2.gls <- function(object, data = NULL, cluster, ...){

    ## ** data
    if(is.null(data)){
        data <- extractData(object, design.matrix = FALSE, as.data.frame = TRUE,
                            envir = parent.env(environment()))
    }

    ## ** endogenous variable
    formula.object <- .getFormula2(object)
    name.Y <- all.vars(stats::update(formula.object, ".~1"))

    ## ** extractors   
    res.cluster <- .getCluster2(object, data = data, cluster = cluster)
    res.param <- .coef2(object)
    res.index <- .getIndexOmega2(object,
                                 param = res.param,
                                 attr.param = attributes(res.param),
                                 name.Y = name.Y,
                                 cluster = res.cluster$cluster,
                                 levels.cluster = res.cluster$levels.cluster,
                                 data = data)
    res.Omega <- .getVarCov2(object,
                             param = res.param,
                             attr.param = attributes(res.param),
                             name.endogenous = res.index$name.endogenous,
                             n.endogenous = res.index$n.endogenous,
                             ref.group = res.index$ref.group)

    ## ** export
    return(c(res.cluster,
             list(param = res.param),
             res.index,
             list(Omega = res.Omega))
           )
}

## * getVarCov2.lme
#' @rdname getVarCov2
#' @export
getVarCov2.lme <- getVarCov2.gls

## * getVarCov2.lvmfit
#' @rdname getVarCov2
#' @export
getVarCov2.lvmfit <- function(object, data = NULL, param = NULL, ...){

    if(inherits(object, "lvmfit2")){
        return(object$sCorrect$Omega)
    }else{
        name.latent <- latent(object)
        n.latent <- length(name.latent)

        ## ** Prepare
        if(is.null(object$conditionalMoment)){
            name.endogenous <- endogenous(object)
            if(is.null(param)){
                param <- coef(object)
            }
            if(is.null(data)){
                data <- as.data.frame(object$data$model.frame)
            }

            object$conditionalMoment <- conditionalMoment(object,
                                                          data = data,
                                                          param = param,
                                                          name.endogenous = name.endogenous,
                                                          name.latent = name.latent,
                                                          first.order = FALSE,
                                                          second.order = FALSE,
                                                          usefit = TRUE)
        }

        ## ** Compute Omega
        if(n.latent>0){
            Omega <- object$conditionalMoment$value$tLambda.tiIB.Psi.iIB %*% object$conditionalMoment$value$Lambda + object$conditionalMoment$value$Sigma
        }else{
            Omega <- object$conditionalMoment$value$Sigma
        }

        return(Omega)
    }
}



##----------------------------------------------------------------------
### getVarCov2.R ends here
