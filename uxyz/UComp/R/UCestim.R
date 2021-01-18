#' @title UCestim
#' @description Estimates and forecasts UC models
#'
#' @details \code{UCestim} estimates and forecasts a time series using an
#' UC model.
#' Standard methods applicable to UComp objects are print, summary, plot,
#' fitted, residuals, logLik, AIC, BIC, coef, predict, tsdiag.
#'
#' @param sys an object of type \code{UComp} created with \code{UC}
#' 
#' @return The same input object with the appropriate fields 
#' filled in, in particular:
#' \item{p}{Estimated parameters}
#' \item{v}{Estimated innovations (white noise in correctly specified models)}
#' \item{yFor}{Forecasted values of output}
#' \item{yForV}{Variance of forecasted values of output}
#' \item{criteria}{Value of criteria for estimated model}
#' 
#' @author Diego J. Pedregal
#' 
#' @seealso \code{\link{UC}}, \code{\link{UCmodel}}, \code{\link{UCvalidate}}, \code{\link{UCfilter}}, 
#'          \code{\link{UCsmooth}}, \code{\link{UCdisturb}}, \code{\link{UCcomponents}},
#'          \code{\link{UChp}}
#'          
#' @examples
#' m1 <- UCsetup(log(AirPassengers))
#' m1 <- UCestim(m1)
#' @rdname UCestim
#' @export
UCestim = function(sys){
    # Clear variables to make new estimation
    sys$table = NA
    sys$hidden$constPar = NA
    # Estimation
    rubbish = c(sys$hidden$d_t, sys$hidden$innVariance, sys$hidden$objFunValue, sys$cLlik, sys$outlier, sys$arma)
    rubbish2 = cbind(sys$hidden$grad, sys$hidden$constPar, sys$hidden$typePar)
    rubbish3 = cbind(sys$hidden$ns, sys$hidden$nPar)
    if (is.ts(sys$u)){
        u = as.numeric(sys$u)
    } else {
        u = sys$u
    }
    nu = dim(u)[2]
    kInitial = dim(u)[1]
    if (nu == 2){
        nu = length(sys$y) + sys$h
        kInitial = 0
    }
    if (is.ts(sys$y)){
        output = UCompC("estimate", as.numeric(sys$y), u, sys$model, sys$periods, sys$rhos,
                    sys$h, sys$tTest, sys$criterion, sys$p, rubbish2, rubbish, sys$verbose, 
                    sys$stepwise, sys$hidden$estimOk, sys$p0, sys$v, sys$yFitV,
                    sys$hidden$nonStationaryTerms, rubbish3, sys$hidden$harmonics,
                    as.vector(sys$criteria), sys$hidden$cycleLimits, 
                    cbind(sys$hidden$beta, sys$hidden$betaV), sys$hidden$typeOutliers)
        fY = frequency(sys$y)
        sY = start(sys$y, frequency = fY)
        # sys$v = ts(output$v, sY, frequency = fY)
        # sys$vV = ts(output$vV, sY, frequency = fY)
        aux = ts(matrix(NA, length(sys$y) + 1, 1), sY, frequency = fY)
        if (length(output$yFor > 0)){
            sys$yFor = ts(output$yFor, end(aux), frequency = fY)
            sys$yForV = ts(output$yForV, end(aux), frequency = fY)
        }
    } else {
        output = UCompC("estimate", sys$y, u, sys$model, sys$periods, sys$rhos,
                        sys$h, sys$tTest, sys$criterion, sys$p, rubbish2, rubbish, sys$verbose, 
                        sys$stepwise, sys$hidden$estimOk, sys$p0, sys$v, sys$yFitV,
                        sys$hidden$nonStationaryTerms, rubbish3, sys$hidden$harmonics,
                        as.vector(sys$criteria), sys$hidden$cycleLimits,
                        cbind(sys$hidden$beta, sys$hidden$betaV), sys$hidden$typeOutliers)
        # sys$v = output$v
        # sys$vV = output$vV
        if (length(output$yFor > 0)){
            sys$yFor = output$yFor
            sys$yForV = output$yForV
        }
    }
    # Converte to R list
    sys$p = output$p
    sys$p0 = output$p0
    # sys$stdP = output$stdP
    if (grepl("?", sys$model, fixed = TRUE)){
        sys$model = output$model
    }
    # sys$hidden$grad = output$grad
    # sys$hidden$constPar = output$constPar
    n = length(sys$p)
    rubbish2 = matrix(output$rubbish2, n, 3)
    sys$hidden$grad = rubbish2[, 1]
    sys$hidden$constPar = rubbish2[, 2]
    sys$hidden$typePar = rubbish2[, 3]
    sys$hidden$cycleLimits = matrix(output$cycleLimits, 
                                    length(output$cycleLimits) / 2, 2)
    sys$hidden$d_t = output$rubbish[1]
    sys$hidden$innVariance = output$rubbish[2]
    sys$hidden$objFunValue = output$rubbish[3]
    betas = matrix(output$betas, length(output$betas) / 2, 2)
    sys$hidden$beta = betas[, 1]
    sys$hidden$betaV = betas[, 2]
    sys$periods = output$periods
    sys$rhos = output$rhos
    sys$hidden$estimOk = output$estimOk
    sys$hidden$nonStationaryTerms = output$nonStationaryTerms
    # sys$hidden$ns = output$ns
    # sys$hidden$nPar = output$nPar
    rubbish3 = matrix(output$rubbish3, 5, 2)
    sys$hidden$ns = rubbish3[, 1]
    sys$hidden$nPar = rubbish3[, 2]
    sys$hidden$harmonics = output$harmonics
    criteria = output$criteria;
    sys$criteria = matrix(criteria, 1, 4)
    colnames(sys$criteria) = c("LLIK", "AIC", "BIC", "AICc")
    if (!is.na(sys$outlier)){
        k = length(output$u) / nu
        nOut = k - kInitial
        if (nOut > 0){
            sys$u = matrix(output$u, k, nu)
            sys$hidden$typeOutliers = matrix(NA, nOut, 2)
            sys$hidden$typeOutliers[, 1] = c(output$typeOutliers);
            for (i in 1 : nOut){
                sys$hidden$typeOutliers[i, 2] = min(which(sys$u[kInitial + i, ] == 1))
            }
        }
    }
    return(sys)
}
    