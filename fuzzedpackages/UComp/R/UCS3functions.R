#' @title print.UComp
#' @description Prints an UComp object
#'
#' @details See help of \code{UC}.
#'
#' @param x Object of class \dQuote{UComp}.
#' @param ... Additional inputs to handle the way to print output.
#' 
#' @author Diego J. Pedregal
#' 
#' @seealso \code{\link{UC}}, \code{\link{UCmodel}}, \code{\link{UCvalidate}}, \code{\link{UCfilter}}, \code{\link{UCsmooth}}, 
#'          \code{\link{UCdisturb}}, \code{\link{UCcomponents}}
#'          
#' @examples
#' y <- log(AirPassengers)
#' m1 <- UCmodel(y, model = "llt/equal/arma(0,0)")
#' print(m1)
#' @noRd
#' @export 
print.UComp = function(x, ...){
    x = UCvalidate(x, TRUE)
}
#' @title summary.UComp
#' @description Prints an UComp object on screen
#'
#' @param object Object of class \dQuote{UComp}.
#' @param ... Additional inputs to function.
#' 
#' @details See help of \code{UC}.
#'
#' @author Diego J. Pedregal
#' 
#' @seealso \code{\link{UC}}, \code{\link{UCmodel}}, \code{\link{UCvalidate}}, \code{\link{UCfilter}}, \code{\link{UCsmooth}}, 
#'          \code{\link{UCdisturb}}, \code{\link{UCcomponents}}
#'          
#' @examples
#' y <- log(AirPassengers)
#' m1 <- UCmodel(y, model = "llt/equal/arma(0,0)")
#' summary(m1)
#' @noRd
#' @export 
summary.UComp = function(object, ...){
    print(object)
}
#' @title plot.UComp
#' @description Plot components of UComp object
#'
#' @details See help of \code{UC}.
#'
#' @param x Object of class \dQuote{UComp}.
#' @param ... Additional inputs to function.
#' 
#' @author Diego J. Pedregal
#' 
#' @seealso \code{\link{UC}}, \code{\link{UCmodel}}, \code{\link{UCvalidate}}, \code{\link{UCfilter}}, \code{\link{UCsmooth}}, 
#'          \code{\link{UCdisturb}}, \code{\link{UCcomponents}}
#'          
#' @examples
#' y <- log(AirPassengers)
#' m1 <- UCmodel(y, model = "llt/equal/arma(0,0)")
#' plot(m1)
#' @noRd
#' @export 
plot.UComp = function(x, ...){
    if (length(x$comp) < 2){
        x = UCcomponents(x)
    }
    if (is.ts(x$comp)){
        plot(x$comp, main = "Time Series Decomposition")
    } else {
        plot(ts(x$comp, frequency = x$periods[1]),
             main = "Time Series Decomposition")
    }
}
#' @title fitted.UComp
#' @description Fitted output values of UComp object
#'
#' @details See help of \code{UC}.
#'
#' @param object Object of class \dQuote{UComp}.
#' @param ... Additional inputs to function.
#' 
#' @author Diego J. Pedregal
#' 
#' @seealso \code{\link{UC}}, \code{\link{UCmodel}}, \code{\link{UCvalidate}}, \code{\link{UCfilter}}, \code{\link{UCsmooth}}, 
#'          \code{\link{UCdisturb}}, \code{\link{UCcomponents}}
#'          
#' @examples
#' y <- log(AirPassengers)
#' m1 <- UCmodel(y, model = "llt/equal/arma(0,0)")
#' fitted(m1)
#' @noRd
#' @export 
fitted.UComp = function(object, ...){
    if (length(object$yFit) < 2){
        object = UCsmooth(object)
    }
    return(object$yFit)
}
#' @title residuals.UComp
#' @description Standardised innovations values of UComp object
#'
#' @details See help of \code{UC}.
#'
#' @param object Object of class \dQuote{UComp}.
#' @param ... Additional inputs to function.
#' 
#' @author Diego J. Pedregal
#' 
#' @seealso \code{\link{UC}}, \code{\link{UCmodel}}, \code{\link{UCvalidate}}, \code{\link{UCfilter}}, \code{\link{UCsmooth}}, 
#'          \code{\link{UCdisturb}}, \code{\link{UCcomponents}}
#'          
#' @examples
#' y <- log(AirPassengers)
#' m1 <- UCmodel(y, model = "llt/equal/arma(0,0)")
#' residuals(m1)
#' @noRd
#' @export 
residuals.UComp = function(object, ...){
    if (length(object$yFit) < 2){
        object = UCfilter(object)
    }
    return(object$v)
}
#' @title logLik.UComp
#' @description Extract log Likelihood value of UComp object
#'
#' @details See help of \code{UC}.
#'
#' @param object Object of class \dQuote{UComp}.
#' @param ... Additional inputs to function.
#' 
#' @author Diego J. Pedregal
#' 
#' @seealso \code{\link{UC}}, \code{\link{UCmodel}}, \code{\link{UCvalidate}}, \code{\link{UCfilter}}, \code{\link{UCsmooth}}, 
#'          \code{\link{UCdisturb}}, \code{\link{UCcomponents}}
#'          
#' @examples
#' y <- log(AirPassengers)
#' m1 <- UCmodel(y, model = "llt/equal/arma(0,0)")
#' logLik(m1)
#' @noRd
#' @export 
logLik.UComp = function(object, ...){
    out = object$criteria[1]
    class(out) = "logLik"
    attr(out, "df") = length(object$p) - 1
    return(out)
}
#' @title AIC.UComp
#' @description Extract AIC value of UComp object
#'
#' @details See help of \code{UC}.
#'
#' @param object Object of class \dQuote{UComp}.
#' @param ... Additional inputs to function.
#' 
#' @author Diego J. Pedregal
#' 
#' @seealso \code{\link{UC}}, \code{\link{UCmodel}}, \code{\link{UCvalidate}}, \code{\link{UCfilter}}, \code{\link{UCsmooth}}, 
#'          \code{\link{UCdisturb}}, \code{\link{UCcomponents}}
#'          
#' @examples
#' y <- log(AirPassengers)
#' m1 <- UCmodel(y, model = "llt/equal/arma(0,0)")
#' logLik(m1)
#' @noRd
#' @export 
AIC.UComp = function(object, ..., k = 2){
    return(object$criteria[2])
}
#' @title BIC.UComp
#' @description Extract BIC (or SBC) value of UComp object
#'
#' @details See help of \code{UC}.
#'
#' @param object Object of class \dQuote{UComp}.
#' @param ... Additional inputs to function.
#' 
#' @author Diego J. Pedregal
#' 
#' @seealso \code{\link{UC}}, \code{\link{UCmodel}}, \code{\link{UCvalidate}}, \code{\link{UCfilter}}, \code{\link{UCsmooth}}, 
#'          \code{\link{UCdisturb}}, \code{\link{UCcomponents}}
#'          
#' @examples
#' y <- log(AirPassengers)
#' m1 <- UCmodel(y, model = "llt/equal/arma(0,0)")
#' logLik(m1)
#' @noRd
#' @export 
BIC.UComp = function(object, ...){
    return(object$criteria[3])
}
#' @title coef.UComp
#' @description Extract model coefficients of UComp object
#'
#' @details See help of \code{UC}.
#'
#' @param object Object of class \dQuote{UComp}.
#' @param ... Additional inputs to function.
#' 
#' @author Diego J. Pedregal
#' 
#' @seealso \code{\link{UC}}, \code{\link{UCmodel}}, \code{\link{UCvalidate}}, \code{\link{UCfilter}}, \code{\link{UCsmooth}}, 
#'          \code{\link{UCdisturb}}, \code{\link{UCcomponents}}
#'          
#' @examples
#' y <- log(AirPassengers)
#' m1 <- UCmodel(y, model = "llt/equal/arma(0,0)")
#' coef(m1)
#' @noRd
#' @export 
coef.UComp = function(object, ...){
    if (length(object$table) < 2){
        object = UCvalidate(object, FALSE)
    }
    # Number of u's
    nu = dim(object$u)[1]
    if (dim(object$u)[2] == 2){
        nu = 0
    }
    parameters = matrix(NA, length(object$p) + nu, 1)
    parametersNames = matrix("", length(parameters))
    # Looking for parameters in output table
    hyphen = 1
    rowM = 1
    i = 1
    while (hyphen < 4){
        rowM = rowM + 1
        if (hyphen > 2){
            lineI = object$table[rowM]
            if (substr(object$table[rowM], 1, 1) != "-"){
                # Parameter name
                namePar = substr(lineI, 1, gregexpr(pattern =':', lineI))
                colonN = nchar(namePar) + 1
                namePar = gsub(" ", "", substr(namePar, 1, colonN - 2), fixed = TRUE)
                # Parameter value
                blanks = gregexpr(pattern =' ', substr(lineI, colonN, nchar(lineI)))
                aux = which(diff(blanks[[1]]) > 1)[1] + 1
                aux2 = gsub("*", "", substr(lineI, colonN, colonN + blanks[[1]][aux]), fixed = TRUE)
                parameters[i, 1] = as.numeric(gsub(" ", "", aux2))
                parametersNames[i] = namePar
                i = i + 1
            } 
            # else {
            #     hyphen = hyphen + 1
            # }
        }
        if (substr(object$table[rowM], 1, 1) == "-"){
            hyphen = hyphen + 1
        }
    }
    rownames(parameters) = parametersNames
    return(parameters)
}
#' @title predict.UComp
#' @description Forecasting using structural Unobseved Components models
#'
#' @details See help of \code{UC}.
#'
#' @param object Object of class \dQuote{UComp}.
#' @param ... Additional inputs to function.
#' 
#' @author Diego J. Pedregal
#' 
#' @return A list with components \code{pred} for the predictions and 
#'         \code{se} for standard errors (if \code{se.fit = TRUE})
#' 
#' @seealso \code{\link{UC}}, \code{\link{UCmodel}}, \code{\link{UCvalidate}}, \code{\link{UCfilter}}, \code{\link{UCsmooth}}, 
#'          \code{\link{UCdisturb}}, \code{\link{UCcomponents}}
#'          
#' @examples
#' y <- log(AirPassengers)
#' m1 <- UCmodel(y, model = "llt/eq/arma(0,0)")
#' f1 <- predict(m1)
#' @noRd
#' @export 
predict.UComp = function(object, ...){
    out = list(pred = object$yFor,
               se = sqrt(object$yForV))
    return(out)
}
#' @title tsdiag.UComp
#' @description Diagnostic plots for UComp objects
#'
#' @details See help of \code{UC}.
#'
#' @param object Object of class \dQuote{UComp}.
#' @param gof.lag Maximum number of lags for pormanteau Ljung-Box test
#' @param ... Additional inputs to function.
#' 
#' @author Diego J. Pedregal
#' 
#' @seealso \code{\link{UC}}, \code{\link{UCmodel}}, \code{\link{UCvalidate}}, \code{\link{UCfilter}}, \code{\link{UCsmooth}}, 
#'          \code{\link{UCdisturb}}, \code{\link{UCcomponents}}
#'          
#' @examples
#' y <- log(AirPassengers)
#' m1 <- UCmodel(y, model = "llt/eq/arma(0,0)")
#' tsdiag(m1)
#' @noRd
#' @export 
tsdiag.UComp = function(object, gof.lag = NULL, ...){
    if (length(object$v) < 2){
        object = UCfilter(object)
    }
    # Number of u's
    nu = dim(object$u)[1]
    if (dim(object$u)[2] == 2){
        nu = 0
    }
    # Number of parameters
    nPar = length(object$p) - 1 + nu
    aux = list(residuals = object$v,
               sigma2 = 1,
               nobs = length(object$y) - nPar,
               coef = object$p,
               x = object$y,
               fitted = object$yFit)
    if (is.null(gof.lag)){
        tsdiag(structure(aux, class = "Arima"))
    } else {
        tsdiag(structure(aux, class = "Arima"), gof.lag)
    }
}

