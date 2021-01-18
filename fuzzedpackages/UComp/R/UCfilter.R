#' @title filter_
#' @description Auxiliar function of \code{UComp} package.
#'
#' @param sys reserved input
#' @param command reserved input
#' 
#' @author Diego J. Pedregal
#' 
#' @noRd
filter_ = function(sys, command){
    if (is.ts(sys$y)){
        y = as.numeric(sys$y)
    } else {
        y = sys$y
    }
    if (is.ts(sys$u)){
        u = as.numeric(sys$u)
    } else {
        u = sys$u
    }
    rubbish = c(sys$hidden$d_t, sys$hidden$innVariance, sys$hidden$objFunValue, sys$cLlik, sys$outlier, sys$arma)
    rubbish2 = cbind(sys$hidden$grad, sys$hidden$constPar, sys$hidden$typePar)
    rubbish3 = cbind(sys$hidden$ns, sys$hidden$nPar)
    output = UCompC(command, y, u, sys$model, sys$periods, sys$rhos,
                    sys$h, sys$tTest, sys$criterion, sys$p, rubbish2, rubbish, sys$verbose, 
                    sys$stepwise, sys$hidden$estimOk, sys$p0, sys$v, sys$yFitV,
                    sys$hidden$nonStationaryTerms, rubbish3, sys$hidden$harmonics,
                    as.vector(sys$criteria), sys$hidden$cycleLimits, 
                    cbind(sys$hidden$beta, sys$hidden$betaV), sys$hidden$typeOutliers)
    # Convert to R list
    # if (command == "disturb"){
    #     sys$eta = output$eta
    #     sys$eps = output$eps
    # } else {
    #     sys$a = output$a
    #     sys$P = output$P
    #     sys$v = output$v
    #     sys$yFitV = output$yFitV
    #     sys$yFit = output$yFit
    # }
    # Re-building matrices to their original sizes
    n = length(sys$y) + sys$h
    # m = length(output$a) / n
    m = sum(sys$hidden$ns)
    if (is.ts(sys$y)){
        fY = frequency(sys$y)
        sY = start(sys$y, frequency = fY)
        # aux = ts(matrix(NA, length(sys$y) + 1, 1), sY, frequency = fY)
        # sys$yFor = ts(output$yFor, end(aux), frequency = fY)
        # sys$yForV = ts(output$yForV, end(aux), frequency = fY)
        if (command == "disturb"){
            n = length(sys$y)
            mEta = length(output$eta) / n
            sys$eta = ts(t(matrix(output$eta, mEta, n)), sY, frequency = fY)
            sys$eps = ts(t(matrix(output$eps, 1, n)), sY, frequency = fY)
        } else {
            sys$a = ts(t(matrix(output$a, m, n)), sY, frequency = fY)
            # sys$aFor = ts(t(matrix(sys$aFor, m, sys$h)), end(aux), frequency = fY)
            sys$P = ts(t(matrix(output$P, m, n)), sY, frequency = fY)
            # sys$PFor = ts(t(matrix(sys$PFor, m, sys$h)), end(aux), frequency = fY)
            sys$yFit = ts(output$yFit, sY, frequency = fY)
            sys$yFitV = ts(output$yFitV, sY, frequency = fY)
            aux = ts(matrix(NA, length(sys$y) - length(output$v) + 1, 1), sY, frequency = fY)
            sys$v = ts(output$v, end(aux), frequency = fY)
        }
    } else {
        if (command == "disturb"){
            n = length(sys$y)
            sys$eta = t(matrix(output$eta, m, n))
            sys$eps = t(matrix(output$eps, 1, n))
        } else {
            sys$a = t(matrix(output$a, m, n))
            # sys$aFor = t(matrix(sys$aFor, m, sys$h))
            sys$P = t(matrix(output$P, m, n))
            # sys$PFor = t(matrix(sys$PFor, m, sys$h))
            sys$yFit = output$yFit
            sys$yFitV = output$yFitV
            sys$v = output$v
        }
    }
    return(sys)
}

#' @title UCfilter
#' @description Runs the Kalman Filter for UC models
#' Standard methods applicable to UComp objects are print, summary, plot,
#' fitted, residuals, logLik, AIC, BIC, coef, predict, tsdiag.
#'
#' @param sys an object of type \code{UComp} created with \code{UC}
#' 
#' @return The same input object with the appropriate fields 
#' filled in, in particular:
#' \item{yFit}{Fitted values of output}
#' \item{yFitV}{Variance of fitted values of output}
#' \item{a}{State estimates}
#' \item{P}{Variance of state estimates}
#' 
#' @author Diego J. Pedregal
#' 
#' @seealso \code{\link{UC}}, \code{\link{UCmodel}}, \code{\link{UCvalidate}}, 
#'          \code{\link{UCsmooth}}, \code{\link{UCdisturb}}, \code{\link{UCcomponents}},
#'          \code{\link{UChp}}
#'          
#' @examples
#' m1 <- UC(log(AirPassengers))
#' m1 <- UCfilter(m1)
#' @rdname UCfilter
#' @export
UCfilter = function(sys){
    return(filter_(sys, "filter"))
}

#' @title UCsmooth
#' @description Runs the Fixed Interval Smoother for UC models
#' Standard methods applicable to UComp objects are print, summary, plot,
#' fitted, residuals, logLik, AIC, BIC, coef, predict, tsdiag.
#'
#' @param sys an object of type \code{UComp} created with \code{UC}
#' 
#' @return The same input object with the appropriate fields 
#' filled in, in particular:
#' \item{yFit}{Fitted values of output}
#' \item{yFitV}{Variance of fitted values of output}
#' \item{a}{State estimates}
#' \item{P}{Variance of state estimates}
#' 
#' @author Diego J. Pedregal
#' 
#' @seealso \code{\link{UC}}, \code{\link{UCmodel}}, \code{\link{UCvalidate}}, \code{\link{UCfilter}}, 
#'          \code{\link{UCdisturb}}, \code{\link{UCcomponents}},
#'          \code{\link{UChp}}
#'          
#' @examples
#' m1 <- UC(log(AirPassengers))
#' m1 <- UCsmooth(m1)
#' @rdname UCsmooth
#' @export
UCsmooth = function(sys){
    return(filter_(sys, "smooth"))
}

#' @title UCdisturb
#' @description Runs the Disturbance Smoother for UC models
#' Standard methods applicable to UComp objects are print, summary, plot,
#' fitted, residuals, logLik, AIC, BIC, coef, predict, tsdiag.
#'
#' @param sys an object of type \code{UComp} created with \code{UC}
#' 
#' @return The same input object with the appropriate fields 
#' filled in, in particular:
#' \item{yFit}{Fitted values of output}
#' \item{yFitV}{Variance of fitted values of output}
#' \item{a}{State estimates}
#' \item{P}{Variance of state estimates}
#' \item{eta}{State perturbations estimates}
#' \item{eps}{Observed perturbations estimates}
#' 
#' @author Diego J. Pedregal
#' 
#' @seealso \code{\link{UC}}, \code{\link{UCmodel}}, \code{\link{UCvalidate}}, \code{\link{UCfilter}}, 
#'          \code{\link{UCsmooth}}, \code{\link{UCcomponents}},
#'          \code{\link{UChp}}
#'          
#' @examples
#' m1 <- UC(log(AirPassengers))
#' m1 <- UCdisturb(m1)
#' @rdname UCdisturb
#' @export
UCdisturb = function(sys){
    return(filter_(sys, "disturb"))
}

#' @title UChp
#' @description Hodrick-Prescott filter estimation
#'
#' @param y A time series object
#' 
#' @param lambda Smoothing constant (default: 1600)
#' 
#' @return The cycle estimation
#' 
#' @author Diego J. Pedregal
#' 
#' @seealso \code{\link{UC}}, \code{\link{UCmodel}}, \code{\link{UCvalidate}}, \code{\link{UCfilter}}, 
#'          \code{\link{UCsmooth}}, \code{\link{UCcomponents}}, \code{\link{UCdisturb}}
#'          
#' @examples
#' cycle <- UChp(USgdp)
#' plot(cycle)
#' @rdname UChp
#' @export
UChp = function(y, lambda = 1600){
    m = UCsetup(y, model = "irw/none/arma(0,0)")
    m$p = c(log(1 / lambda) / 2, 0)
    m = UCcomponents(m)
    return(y - m$comp[, 1])
}
