### residuals2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov  8 2017 (09:05) 
## Version: 
## Last-Updated: feb 11 2019 (13:25) 
##           By: Brice Ozenne
##     Update #: 935
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * documentation - residuals2
#' @title Extract Corrected Residuals
#' @description Extract correct residuals from a gaussian linear model.
#' @name residuals2
#' 
#' @param object a \code{lm2}, \code{gls2}, \code{lme2}, or \code{lvmfit2} object.
#' @param type [character] the type of residual to extract:
#' \code{"response"} for raw residuals,
#' \code{"studentized"} for studentized residuals,
#' \code{"normalized"} for normalized residuals.
#' @param param [named numeric vector] the fitted parameters.
#' @param data [data.frame] the data set.
#'
#' @seealso \code{\link{sCorrect}} to obtain \code{lm2}, \code{gls2}, \code{lme2}, or \code{lvmfit2} objects.
#'
#' @details If argument \code{p} or \code{data} is not null, then the small sample size correction is recomputed to correct the residuals. \cr
#'
#' The raw residuals are defined by  observation minus the fitted value:
#' \deqn{
#' \varepsilon = (Y_1 - \mu_1, ..., Y_m - \mu_m)
#' }
#' The studentized residuals divided the raw residuals relative to each endogenous variable by the modeled variance of the endogenous variable.
#' \deqn{
#' \varepsilon_{stud} =(\frac{Y_1 - \mu_1}{\sigma_1}, ..., \frac{Y_m - \mu_m}{\sigma_m})
#' }
#' The normalized residuals multiply the raw residuals by the inverse of the square root of the modeled residual variance covariance matrix.
#' \deqn{
#' \varepsilon_{norm} = \varepsilon \Omega^{-1/2}
#' }
#' @return a matrix containing the residuals relative to each sample (in rows)
#' and each endogenous variable (in column).
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
#' sCorrect(e.lm) <- TRUE
#' 
#' sigma(e.lm)^2
#' mean(residuals(e.lm)^2)
#' mean(residuals2(e.lm)^2)
#' 
#' ## latent variable model
#' e.lvm <- estimate(m, data = d)
#' sCorrect(e.lvm) <- TRUE
#' mean(residuals2(e.lvm)^2)
#'
#' @concept small sample inference
#' @export
`residuals2` <-
    function(object, param, data, type) UseMethod("residuals2")

## * residuals2.lm2
#' @rdname residuals2
#' @export
residuals2.lm2 <- function(object, param = NULL, data = NULL, type = "response"){

    type <- match.arg(type, choices = c("response","studentized","normalized"), several.ok = FALSE)

    if(!is.null(param) || !is.null(data)){
        args <- object$sCorrect$args
        args$df <- FALSE
        object$sCorrect <- do.call(sCorrect,
                                   args = c(list(object, param = param, data = data),
                                            args))
    }
    if(type=="response"){
        residuals <- object$sCorrect$residuals
    }else if(type=="studentized"){
        residuals <- sweep(object$sCorrect$residuals,
                           STATS = sqrt(diag(object$sCorrect$Omega)),
                           FUN = "/",
                           MARGIN = 2)
        ## object$sCorrect$residuals/residuals
    }else if(type=="normalized"){
        residuals <- object$sCorrect$residuals %*% matrixPower(object$sCorrect$Omega, symmetric = TRUE, power = -1/2)
        colnames(residuals) <- colnames(object$sCorrect$residuals)
        ## object$sCorrect$residuals/residuals
        ## var(residuals)        
    }
    return(residuals)
}

## * residuals2.gls2
#' @rdname residuals2
#' @export
residuals2.gls2 <- residuals2.lm2

## * residuals2.lme2
#' @rdname residuals2
#' @export
residuals2.lme2 <- residuals2.lm2

## * residuals2.lvmfit2
#' @rdname residuals2
#' @export
residuals2.lvmfit2 <- residuals2.lm2


##----------------------------------------------------------------------
### residuals2.R ends here
