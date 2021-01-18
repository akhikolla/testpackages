### score2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 12 2017 (16:43) 
## Version: 
## last-updated: jul 31 2020 (10:44) 
##           By: Brice Ozenne
##     Update #: 2263
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * Documentation - score2
#' @title  Extract the Individual Score
#' @description  Extract the individual score from a Gaussian linear model.
#' @name score2
#'
#' @param object a linear model or a latent variable model
#' @param param [optional] the fitted parameters.
#' @param data [optional] the data set.
#' @param bias.correct [logical] should the standard errors of the coefficients be corrected for small sample bias? Only relevant if the \code{sCorrect} function has not yet be applied to the object.
#' @param ... arguments to be passed to \code{sCorrect}.
#'
#' @details If argument \code{p} or \code{data} is not null, then the small sample size correction is recomputed to correct the influence function.
#'
#' @seealso \code{\link{sCorrect}} to obtain \code{lm2}, \code{gls2}, \code{lme2}, or \code{lvmfit2} objects.
#'
#' @return A matrix containing the score relative to each sample (in rows)
#' and each model coefficient (in columns).
#' 
#' @examples
#' n <- 5e1
#' p <- 3
#' X.name <- paste0("X",1:p)
#' link.lvm <- paste0("Y~",X.name)
#' formula.lvm <- as.formula(paste0("Y~",paste0(X.name,collapse="+")))
#'
#' m <- lvm(formula.lvm)
#' distribution(m,~Id) <- Sequence.lvm(0)
#' set.seed(10)
#' d <- lava::sim(m,n)
#'
#' ## linear model
#' e.lm <- lm(formula.lvm,data=d)
#' score.tempo <- score2(e.lm, bias.correct = FALSE)
#' colMeans(score.tempo)
#'
#' ## latent variable model
#' e.lvm <- estimate(lvm(formula.lvm),data=d)
#' score.tempo <- score2(e.lvm, bias.correct = FALSE)
#' range(score.tempo-score(e.lvm, indiv = TRUE))
#'
#' @concept small sample inference
#' @export
`score2` <-
  function(object, ...) UseMethod("score2")

## * score2.lm
#' @rdname score2
#' @export
score2.lm <- function(object, param = NULL, data = NULL, bias.correct = TRUE, ...){
    sCorrect(object, param = param, data = data, df = FALSE, ...) <- bias.correct

    ### ** export
    return(object$sCorrect$score)
}

## * score2.gls
#' @rdname score2
#' @export
score2.gls <- score2.lm

## * score2.lme
#' @rdname score2
#' @export
score2.lme <- score2.lm

## * score2.lvmfit
#' @rdname score2
#' @export
score2.lvmfit <- score2.lm

## * score2.lm2
#' @rdname score2
#' @export
score2.lm2 <- function(object, param = NULL, data = NULL, ...){

    ### ** compute the score
    if(!is.null(param) || !is.null(data)){
        args <- object$sCorrect$args
        args$df <- FALSE
        object$sCorrect <- do.call(sCorrect,
                                   args = c(list(object, param = param, data = data),
                                            args))
    }

    ### ** export
    return(object$sCorrect$score)

}

## * score2.gls2
#' @rdname score2
#' @export
score2.gls2 <- score2.lm2

## * score2.lme2
#' @rdname score2
#' @export
score2.lme2 <- score2.lm2

## * score2.lvmfit
#' @rdname score2
#' @export
score2.lvmfit2 <- score2.lm2

## * .score2
#' @title Compute the Corrected Score.
#' @description Compute the corrected score when there is no missing value.
#' @name score2-internal
#' 
#' @param n.cluster [integer >0] the number of observations.
#' 
#' @keywords internal
.score2 <- function(epsilon, Omega, OmegaM1, dmu, dOmega,                    
                    name.param, name.meanparam, name.varparam,
                    index.Omega, n.cluster, indiv){

### ** Prepare
    test.global <- is.null(index.Omega)
    out.score <- matrix(0, nrow = n.cluster, ncol = length(name.param),
                        dimnames = list(NULL,name.param))
            
### ** global
    if(test.global){
        epsilon.OmegaM1 <- epsilon %*% OmegaM1

        ## *** Compute score relative to the mean coefficients
        for(iP in name.meanparam){ # iP <- 1
            out.score[,iP] <- out.score[,iP] + rowSums(dmu[[iP]] * epsilon.OmegaM1)
        }
        
        ## *** Compute score relative to the variance-covariance coefficients
        for(iP in name.varparam){ # iP <- 1
            term2 <- - 1/2 * tr(OmegaM1 %*% dOmega[[iP]])            
            term3 <- 1/2 * rowSums(epsilon.OmegaM1 %*% dOmega[[iP]] * epsilon.OmegaM1)
            out.score[,iP] <- out.score[,iP] + as.double(term2) + term3
        }        
    }


### ** individual specific
    if(!test.global){

        for(iC in 1:n.cluster){
            iIndex <- index.Omega[[iC]]
            iEpsilon.OmegaM1 <- OmegaM1[[iC]] %*% cbind(epsilon[iC,iIndex])


            ## *** Compute score relative to the mean coefficients
            for(iP in name.meanparam){ # iP <- name.meanparam[1]
                out.score[iC,iP] <- out.score[iC,iP] + dmu[[iP]][iC,iIndex] %*% iEpsilon.OmegaM1
            }

            ## *** Compute score relative to the variance-covariance coefficients
            for(iP in name.varparam){ # iP <- name.varparam[1]
                term2 <- - 1/2 * tr(OmegaM1[[iC]] %*% dOmega[[iP]][iIndex,iIndex,drop=FALSE])
                term3 <- 1/2 * sum(iEpsilon.OmegaM1 * dOmega[[iP]][iIndex,iIndex,drop=FALSE] %*% iEpsilon.OmegaM1)
                out.score[iC,iP] <- out.score[iC,iP] + as.double(term2) + term3 
            }
        }
        
    }

    ### ** export
    if(indiv==FALSE){
        out.score <- colSums(out.score)
    }
    return(out.score)
}


#----------------------------------------------------------------------
### score2.R ends her
