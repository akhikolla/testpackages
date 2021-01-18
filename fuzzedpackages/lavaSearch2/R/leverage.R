### leverage.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 19 2018 (17:58) 
## Version: 
## Last-Updated: feb 11 2019 (13:26) 
##           By: Brice Ozenne
##     Update #: 36
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * documentation - leverage2
#' @title Extract Leverage Values
#' @description Extract leverage values from a Gaussian linear model. 
#' @name leverage2
#' 
#' @param object a \code{lm2}, \code{gls2}, \code{lme2}, or \code{lvmfit2} object.
#' @param param [optional] the fitted parameters.
#' @param data [optional] the data set.
#' @param ... arguments to be passed to \code{sCorrect}.
#'
#' @details The leverage are defined as the partial derivative of the fitted values with respect to the observations.
#' \deqn{
#' leverage_i = \frac{\partial \hat{Y}_i}{\partial Y_i}
#' }
#' See Wei et al. (1998). \cr \cr
#' 
#' If argument \code{p} or \code{data} is not null, then the small sample size correction is recomputed to correct the residuals.
#'
#' @seealso \code{\link{sCorrect}} to obtain \code{lm2}, \code{gls2}, \code{lme2}, or \code{lvmfit2} objects.
#' 
#' @return a matrix containing the leverage relative to each sample (in rows)
#' and each endogenous variable (in column).
#'
#' @references Bo-Cheng Wei et al., Generalized Leverage and its applications (1998), Scandinavian Journal of Statistics 25:1:25-37.
#' 
#' @examples
#' ## simulate data
#' set.seed(10)
#' m <- lvm(Y1~eta,Y2~eta,Y3~eta)
#' latent(m) <- ~eta
#' d <- lava::sim(m,20, latent = FALSE)
#'
#' ## standard linear model
#' e.lm <- lm(Y1~Y2, data = d)
#'
#' sCorrect(e.lm) <- TRUE
#' range(as.double(leverage2(e.lm)) - influence(e.lm)$hat)
#'
#' ## latent variable model
#' e.lvm <- estimate(m, data = d)
#' sCorrect(e.lvm) <- TRUE
#' leverage2(e.lvm)
#' 
#' @concept small sample inference
#' @export
`leverage2` <-
    function(object, ...) UseMethod("leverage2")

## * leverage2.lm
#' @rdname leverage2
#' @export
leverage2.lm <- function(object, param = NULL, data = NULL, ...){

    sCorrect(object, param = param, data = data, df = FALSE, ...) <- TRUE

    ### ** export
    return(object$sCorrect$leverage)
}

## * leverage2.gls
#' @rdname leverage2
#' @export
leverage2.gls <- leverage2.lm

## * leverage2.lme
#' @rdname leverage2
#' @export
leverage2.lme <- leverage2.lm

## * leverage2.lvmfit
#' @rdname leverage2
#' @export
leverage2.lvmfit <- leverage2.lm

## * leverage2.lm2
#' @rdname leverage2
#' @export
leverage2.lm2 <- function(object, param = NULL, data = NULL, ...){

    if(!is.null(param) || !is.null(data)){
        args <- object$sCorrect$args
        args$df <- FALSE
        object$sCorrect <- do.call(sCorrect,
                                   args = c(list(object, param = param, data = data),
                                            args))
    }
    
    return(object$sCorrect$leverage)   
}

## * leverage2.gls2
#' @rdname leverage2
#' @export
leverage2.gls2 <- leverage2.lm2

## * leverage2.lme2
#' @rdname leverage2
#' @export
leverage2.lme2 <- leverage2.lm2

## * leverage2.lvmfit2
#' @rdname leverage2
#' @export
leverage2.lvmfit2 <- leverage2.lm2



##----------------------------------------------------------------------
### leverage.R ends here
