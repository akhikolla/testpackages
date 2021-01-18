#' @title UCcomponents
#' @description Estimates components of UC models
#' Standard methods applicable to UComp objects are print, summary, plot,
#' fitted, residuals, logLik, AIC, BIC, coef, predict, tsdiag.
#'
#' @param sys an object of type \code{UComp} created with \code{UC}
#' 
#' @return The same input object with the appropriate fields 
#' filled in, in particular:
#' \item{comp}{Estimated components in matrix form}
#' \item{compV}{Estimated components variance in matrix form}
#' 
#' @author Diego J. Pedregal
#' 
#' @seealso \code{\link{UC}}, \code{\link{UCmodel}}, \code{\link{UCvalidate}}, \code{\link{UCfilter}}, 
#'          \code{\link{UCsmooth}}, \code{\link{UCdisturb}},
#'          \code{\link{UChp}}
#'          
#' @examples
#' m1 <- UC(log(AirPassengers))
#' m1 <- UCcomponents(m1)
#' @rdname UCcomponents
#' @export
UCcomponents= function(sys){
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
    output = UCompC("components", y, u, sys$model, sys$periods, sys$rhos,
                    sys$h, sys$tTest, sys$criterion, sys$p, rubbish2, rubbish, sys$verbose, 
                    sys$stepwise, sys$hidden$estimOk, sys$p0, sys$v, sys$yFitV,
                    sys$hidden$nonStationaryTerms, rubbish3, sys$hidden$harmonics,
                    as.vector(sys$criteria), sys$hidden$cycleLimits, 
                    cbind(sys$hidden$beta, sys$hidden$betaV), sys$hidden$typeOutliers)
    # Convert to R list
    sys$comp = output$comp
    sys$compV = output$compV
    m = output$m  # + nCycles
    if (dim(u)[1] == 1 && dim(u)[2] == 2){
        k = 0
    } else {
        k = dim(u)[1]
    }
    nCycles = m - k - 4
    # Re-building matrices to their original sizes
    n = length(sys$comp) / m
    if (is.ts(sys$y)){
        sys$comp = ts(t(matrix(sys$comp, m, n)), start(sys$y, frequency = frequency(sys$y)), frequency = frequency(sys$y))
        sys$compV = ts(t(matrix(sys$compV, m, n)), start(sys$y, frequency = frequency(sys$y)), frequency = frequency(sys$y))
    } else {
        sys$comp = t(matrix(sys$comp, m, n))
        sys$compV = t(matrix(sys$compV, m, n))
    }
    namesComp = c("Trend", "Slope", "Seasonal", "Irregular")
    if (nCycles > 0){
        for (i in 1 : nCycles){
            namesComp = c(namesComp, paste0("Cycle", i))
        }
    }
    # Inputs names
    if (k > 0){
        nOut = 0;
        if (sys$hidden$typeOutliers[1, 2] != -1){
            nOut = dim(sys$hidden$typeOutliers)[1]
        }
        nU = k - nOut
        if (nU > 0){
            for (i in 1 : nU){
                namesComp = c(namesComp, paste0("Exogenous", i))
            }
        }
        if (nOut > 0){
            for (i in 1 : nOut){
                namei = "AO"
                if (sys$hidden$typeOutliers[i, 1] == 1){
                    namei = "LS"
                } else if (sys$hidden$typeOutliers[i, 1] == 2){
                    namei = "SC"
                }
                namesComp = c(namesComp, paste0(namei, sys$hidden$typeOutliers[i, 2]))
            }
        }
    }
    colnames(sys$comp) = namesComp
    colnames(sys$compV) = namesComp
    return(sys)
}
    