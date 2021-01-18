### vcov2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 12 2018 (16:38) 
## Version: 
## Last-Updated: jul 31 2020 (10:44) 
##           By: Brice Ozenne
##     Update #: 13
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation - vcov2
#' @title  Extract the Variance Covariance Matrix of the Model Parameters
#' @description  Extract the variance covariance matrix of the model parameters from a Gaussian linear model.
#' @name vcov2
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
#' @return A matrix.
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
#' vcov.tempo <- vcov2(e.lm, bias.correct = TRUE)
#' vcov.tempo[rownames(vcov(e.lm)),colnames(vcov(e.lm))]/vcov(e.lm)
#'
#' ## latent variable model
#' e.lvm <- estimate(lvm(formula.lvm),data=d)
#' vcov.tempo <- vcov2(e.lvm, bias.correct = FALSE)
#' vcov.tempo/vcov(e.lvm)
#'
#' @concept small sample inference
#' @export
`vcov2` <-
  function(object, ...) UseMethod("vcov2")

## * vcov2.lm
#' @rdname vcov2
#' @export
vcov2.lm <- function(object, param = NULL, data = NULL, bias.correct = TRUE, ...){

    sCorrect(object, param = param, data = data, df = FALSE, ...) <- bias.correct

    ### ** export
    return(object$sCorrect$vcov.param)
}

## * vcov2.gls
#' @rdname vcov2
#' @export
vcov2.gls <- vcov2.lm

## * vcov2.lme
#' @rdname vcov2
#' @export
vcov2.lme <- vcov2.lm

## * vcov2.lvmfit
#' @rdname vcov2
#' @export
vcov2.lvmfit <- vcov2.lm

## * vcov2.lm2
#' @rdname vcov2
#' @export
vcov2.lm2 <- function(object, param = NULL, data = NULL, ...){

    ### ** compute the score
    if(!is.null(param) || !is.null(data)){
        args <- object$sCorrect$args
        args$df <- FALSE
        object$sCorrect <- do.call(sCorrect,
                                   args = c(list(object, param = param, data = data),
                                            args))
    }

    ### ** export
    return(object$sCorrect$vcov.param)

}

## * vcov2.gls
#' @rdname vcov2
#' @export
vcov2.gls2 <- vcov2.lm2

## * vcov2.lme
#' @rdname vcov2
#' @export
vcov2.lme2 <- vcov2.lm2

## * vcov2.lvmfit
#' @rdname vcov2
#' @export
vcov2.lvmfit2 <- vcov2.lm2


##----------------------------------------------------------------------
### vcov2.R ends here
