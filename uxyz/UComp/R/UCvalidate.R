#' @title UCvalidate
#' @description Shows a table of estimation and diagnostics results for UC models.
#' Equivalent to print or summary.
#'
#' @param sys an object of type \code{UComp} created with \code{UC}
#' @param printScreen print to screen or just return output table
#' 
#' @return The same input object with the appropriate fields 
#' filled in, in particular:
#' \item{table}{Estimation and validation table}
#' 
#' @author Diego J. Pedregal
#' 
#' @seealso \code{\link{UC}}, \code{\link{UCmodel}}, \code{\link{UCfilter}}, 
#'          \code{\link{UCsmooth}}, \code{\link{UCdisturb}}, \code{\link{UCcomponents}},
#'          \code{\link{UChp}}
#'          
#' @examples
#' m1 <- UC(log(AirPassengers))
#' m1 <- UCvalidate(m1)
#' @rdname UCvalidate
#' @export
UCvalidate = function(sys, printScreen = TRUE){
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
    if (any(is.na(sys$p))){
        if (printScreen){
            print(sys$table)
        }
        return(sys)
    }
    # Converte to R list
    rubbish = c(sys$hidden$d_t, sys$hidden$innVariance, sys$hidden$objFunValue, sys$cLlik, sys$outlier, sys$arma)
    rubbish2 = cbind(sys$hidden$grad, sys$hidden$constPar, sys$hidden$typePar)
    rubbish3 = cbind(sys$hidden$ns, sys$hidden$nPar)
    output = UCompC("validate", y, u, sys$model, sys$periods, sys$rhos,
                    sys$h, sys$tTest, sys$criterion, sys$p, rubbish2, rubbish, sys$verbose, 
                    sys$stepwise, sys$hidden$estimOk, sys$p0, sys$v, sys$yFitV,
                    sys$hidden$nonStationaryTerms, rubbish3, sys$hidden$harmonics,
                    as.vector(sys$criteria), sys$hidden$cycleLimits, 
                    cbind(sys$hidden$beta, sys$hidden$betaV), sys$hidden$typeOutliers)
    sys$table = output$table
    if (is.ts(sys$y)){
        fY = frequency(sys$y)
        sY = start(sys$y, frequency = fY)
        aux = ts(matrix(NA, length(sys$y) - length(output$v) + 1, 1), sY, frequency = fY)
        sys$v = ts(output$v, end(aux), frequency = fY)
    } else {
        sys$v = output$v
    }
    if (printScreen){
        cat(output$table)
    }
    return(sys)
}
    