#' @title noModel
#' @description Auxiliar function of \dQuote{UComp} library.
#'
#' @param model reserved input
#' @param periods reserved input
#' 
#' @author Diego J. Pedregal
#' 
#' @noRd
noModel = function(model, periods){
    comps = strsplit(tolower(model), "/")
    mT = substr(model, 1, 1)
    mC = substr(comps[[1]][2], 1, 1)
    if (periods[1] < 2){
        mS = "n"
    } else {
        mS = substr(comps[[1]][3], 1, 1)
    }
    mI = substr(comps[[1]][4], 1, 1)
    if (mT == "n" && mC == "n" && mS == "n" && mI == "n")
        return(TRUE)
    else
        return(FALSE)
}
#' @title UCsetup
#' @description Sets up UC general univariate models
#'
#' @details See help of \code{UC}.
#'
#' @param y a time series to forecast (it may be either a numerical vector or
#' a time series object). This is the only input required. If a vector, the additional
#' input \code{periods} should be supplied compulsorily (see below).
#' @param u a matrix of input time series. If 
#' the output wanted to be forecast, matrix \code{u} should contain future values for inputs.
#' @param model the model to estimate. It is a single string indicating the type of 
#' model for each component. It allows two formats "trend/seasonal/irregular" or 
#' "trend/cycle/seasonal/irregular". The possibilities available for each component are:
#' \itemize{
#' \item Trend: ? / none / rw / irw / llt / dt; 
#' 
#' \item Seasonal: ? / none / equal / different;
#' 
#' \item Irregular: ? / none / arma(0, 0) / arma(p, q) - with p and q integer positive orders;
#'     
#' \item Cycles: ? / none / combination of positive or negative numbers. 
#' 
#' Positive numbers fix
#' the period of the cycle while negative values estimate the period taking as initial
#' condition the absolute value of the period supplied.
#' Several cycles with positive or negative values are possible
#' and if a question mark is included, the model test for the existence of the cycles
#' specified (check the examples below).
#' }
#' @param outlier critical level of outlier tests. If NA it does not carry out any 
#' outlier detection (default). A negative value indicates critical minimum t test for one run of
#' outlier detection after identification. A positive value indicates the critical minium
#' t test for outlier detection in any model during identification.
#' @param stepwise stepwise identification procedure (TRUE / FALSE).
#' @param tTest augmented Dickey Fuller test for unit roots (TRUE / FALSE). The number of models
#' to search for is reduced, depending on the result of this test.
#' @param p0 initial condition for parameter estimates.
#' @param h forecast horizon. If the model includes inputs h is not used, the lenght of u is used instead.
#' @param criterion information criterion for identification ("aic", "bic" or "aicc").
#' @param periods vector of fundamental period and harmonics.
#' @param verbose intermediate results shown about progress of estimation (TRUE / FALSE).
#' @param arma check for arma models for irregular components (TRUE / FALSE).
#' @param cLlik reserved input
#' 
#' @author Diego J. Pedregal
#' 
#' @return An object of class \code{UComp}. See \code{UCmodel}.
#' Standard methods applicable to UComp objects are print, summary, plot,
#' fitted, residuals, logLik, AIC, BIC, coef, predict, tsdiag.
#' 
#' @seealso \code{\link{UC}}, \code{\link{UCmodel}}, \code{\link{UCvalidate}}, \code{\link{UCfilter}}, \code{\link{UCsmooth}}, 
#'          \code{\link{UCdisturb}}, \code{\link{UCcomponents}},
#'          \code{\link{UChp}}
#'          
#' @examples
#' y <- log(AirPassengers)
#' m1 <- UCsetup(y)
#' m1 <- UCsetup(y,  model = "llt/equal/arma(0,0)")
#' m1 <- UCsetup(y,  outlier = 4)
#' @rdname UCsetup
#' @export
UCsetup = function(y, u = NULL, model = "?/none/?/?", h = NA, outlier = NA, tTest = FALSE, criterion = "aic",
                   periods = NA, verbose = FALSE, stepwise = FALSE, p0 = NA, cLlik = TRUE, arma = TRUE){
    rhos = NA
    p = NA
    # Converting y vector to matrix
    # if (is.vector(y)){
    #     y = matrix(y, length(y), 1)
    # }
    # Converting u vector to matrix
    if (is.vector(u)){
        u = matrix(u, 1, length(u))
    }
    n = length(y)
    if (is.null(u)){
        u = matrix(0, 1, 2)
    } else {
        if (dim(u)[2] < n){
            u = t(u)
        }
        if (length(y) > dim(u)[2]){
            stop("Length of output data never could be greater than length of inputs!!!")
        }
        h = dim(u)[2] - n
    }
    # Removing nans at beginning or end
    if (is.na(y[1]) || is.na(y[n])){
        ind = which(!is.na(y))
        minInd = min(ind)
        maxInd = max(ind)
        if (is.ts(y)){
            freq = frequency(y)
            starty = start(y, frequency = freq)
            if (minInd > 1){
                aux = ts(matrix(NA, minInd, 1), starty, frequency = freq)
                ini = end(aux)
            } else {
                ini = starty
            }
            if (maxInd < n){
                aux = ts(matrix(NA, maxInd, 1), starty, frequency = freq)
                fin = end(aux)
            } else {
                fin = end(y)
            }
            y = ts(y[minInd : maxInd], start = ini, end = fin, frequency = freq)
        } else {
            y = y[minInd : maxInd]
        }
        if (dim(u)[2] > 2){
            u = u[, minInd : maxInd]
        }
    }
    # Checking periods
    if (is.ts(y) && is.na(periods) && frequency(y) > 1){
        periods = frequency(y) / (1 : floor(frequency(y) / 2))
    } else if (is.ts(y) && is.na(periods)){
        periods = 1
    } else if (!is.ts(y) && is.na(periods)){
        stop("Input \"periods\" should be supplied!!")
    }
    # If period == 1 (annual) then change seasonal model to "none"
    if (periods[1] == 1){
        comps = strsplit(model, "/")
        if (length(comps[[1]]) == 3){
            model = paste(comps[[1]][1], "/none/", comps[[1]][3], sep = "")
        } else {
            model = paste(comps[[1]][1], "/", comps[[1]][2], "/none/", comps[[1]][4], sep = "")
        }
    }
    # Adding cycle in case of T/S/I model specification
    nComp = length(gregexpr("/", model)[[1]])
    if (nComp == 2){
        model = sub("/", "/none/", model)
    }
    # Checking model
    model = tolower(model)
    if (noModel(model, periods))
        stop("No model specified!!!")
    if (any(utf8ToInt(model) == utf8ToInt("?")) && !is.na(p0)){
        p0 = NA;
    }
    if (any(utf8ToInt(model) == utf8ToInt("?")) && !is.na(p)){
        p = NA;
    }
    if (grepl("arma", model) && !grepl('(', model, fixed = TRUE))
        model = paste(model, "(0,0)", sep = "")
    if (grepl("arma", model) && !grepl(')', substr(model, nchar(model) - 1, nchar(model)), fixed = TRUE))
        model = paste(model, ")", sep = "")
    if (nComp < 2 || nComp > 3)
        stop("Incorrect number of components!!!")
    comps = strsplit(model, "/")
    if (regexpr(substr(comps[[1]][1], 1, 1), "?nirld")[1] < 0)
        stop("Trend model not correct!!!")
    if (regexpr(substr(comps[[1]][2], 1, 1), "?n+-0123456789")[1] < 0)
        stop("Cycle model not correct!!!")
    if (regexpr(substr(comps[[1]][3], 1, 1), "?ned")[1] < 0)
        stop("Seasonal model not correct!!!")
    if (regexpr(substr(comps[[1]][4], 1, 1), "?na")[1] < 0)
        stop("Irregular model not correct!!!")
    # Checking horizon
    if (is.na(h)){
        if (is.ts(y)){
            per = frequency(y)
        } else {
            per = max(periods)
        }
        if (per == 1){
            per = 5
        }
        h = 2 * per
    }    
    # Set rhos
    if (is.na(rhos)){
        rhos = matrix(1, length(periods))
    }
    # Checking cycle
    mC0 = strsplit(tolower(model), "/")[[1]][2]
    mC = mC0
    if (mC0 == "?"){
        mC = paste0(toString(-4 * frequency(y)), "?")
    } else if (substr(mC0, 1, 1) != "+" && substr(mC0, 1, 1) != "-" && substr(mC0, 1, 1) != "n"){
        mC = paste0("+", mC0)
    }
    model = sub(paste0("/", mC0, "/"), paste0("/", mC, "/"), model, fixed = TRUE)
    out = list(y = y,
               u = u,
               model = model,
               h = h,
               # Outputs
               comp = NA,
               compV = NA,
               p = p,
               v= NA,
               yFit = NA,
               yFor = NA,
               yFitV = NA,
               yForV = NA,
               a = NA,
               P = NA,
               eta = NA,
               eps = NA,
               table = NA,
               # Other less important
               arma =  arma,
               outlier = -outlier,
               tTest = tTest,
               criterion = criterion,
               periods = periods,
               rhos = rhos,
               verbose = verbose,
               stepwise = stepwise,
               p0 = p0,
               cLlik = cLlik,
               criteria = NA,
               # Other variables
               hidden = list(grad = NA,
                             d_t = 0,
                             estimOk = "Not estimated",
                             objFunValue = 0,
                             innVariance = 1,
                             nonStationaryTerms = NA,
                             ns = NA,
                             nPar = NA,
                             harmonics = 0,
                             constPar = NA,
                             typePar = NA,
                             cycleLimits = NA,
                             typeOutliers = matrix(-1, 1, 2),
                             beta = NA,
                             betaV = NA))
    return(structure(out, class = "UComp"))
}

#' @title UCmodel
#' @description Estimates and forecasts UC general univariate models
#'
#' @details \code{UCmodel} is a function for modelling and forecasting univariate
#' time series according to Unobserved Components models (UC). 
#' It sets up the model with a number of control variables that
#' govern the way the rest of functions in the package will work. It also estimates 
#' the model parameters by Maximum Likelihood and forecasts the data.
#' Standard methods applicable to UComp objects are print, summary, plot,
#' fitted, residuals, logLik, AIC, BIC, coef, predict, tsdiag.
#'
#' @inheritParams UCsetup
#' 
#' @return An object of class \code{UComp}. It is a list with fields including all the inputs and
#'         the fields listed below as outputs. All the functions in this package fill in
#'         part of the fields of any \code{UComp} object as specified in what follows (function 
#'         \code{UC} fills in all of them at once):
#' 
#' After running \code{UCmodel} or \code{UCestim}:
#' \item{p}{Estimated parameters}
#' \item{v}{Estimated innovations (white noise in correctly specified models)}
#' \item{yFor}{Forecasted values of output}
#' \item{yForV}{Variance of forecasted values of output}
#' \item{criteria}{Value of criteria for estimated model}
#' 
#' After running \code{UCvalidate}:
#' \item{table}{Estimation and validation table}
#' 
#' After running \code{UCcomponents}:
#' \item{comp}{Estimated components in matrix form}
#' \item{compV}{Estimated components variance in matrix form}
#' 
#' After running \code{UCfilter}, \code{UCsmooth} or  \code{UCdisturb}:
#' \item{yFit}{Fitted values of output}
#' \item{yFitV}{Variance of fitted values of output}
#' \item{a}{State estimates}
#' \item{P}{Variance of state estimates}
#' \item{aFor}{Forecasts of states}
#' \item{PFor}{Forecasts of states variances}
#' 
#' After running \code{UCdisturb}:
#' \item{eta}{State perturbations estimates}
#' \item{eps}{Observed perturbations estimates}
#' 
#' @author Diego J. Pedregal
#' 
#' @seealso \code{\link{UC}}, \code{\link{UCvalidate}}, \code{\link{UCfilter}}, \code{\link{UCsmooth}}, 
#'          \code{\link{UCdisturb}}, \code{\link{UCcomponents}},
#'          \code{\link{UChp}}
#'          
#' @examples
#' y <- log(AirPassengers)
#' m1 <- UCmodel(y)
#' m1 <- UCmodel(y, model = "llt/equal/arma(0,0)")
#' @rdname UCmodel
#' @export
UCmodel = function(y, u = NULL, model = "?/none/?/?", h = NA, outlier = NA, tTest = FALSE, criterion = "aic",
                   periods = NA, verbose = FALSE, stepwise = FALSE, p0 = NA, cLlik = TRUE, arma = TRUE){
    m1 = UCsetup(y, u, model, h, outlier, tTest, criterion, 
                 periods, verbose, stepwise, p0, cLlik, arma)
    m1 = UCestim(m1)
    return(m1)
}
#' @title UC
#' @description Runs all relevant functions for UC modelling
#'
#' @details See help of \code{UCmodel}.
#'
#' @inheritParams UCsetup
#' 
#' @author Diego J. Pedregal
#' 
#' @return An object of class \code{UComp}. See \code{UC}.
#' Standard methods applicable to UComp objects are print, summary, plot,
#' fitted, residuals, logLik, AIC, BIC, coef, predict, tsdiag.
#' 
#' @seealso \code{\link{UC}}, \code{\link{UCvalidate}}, \code{\link{UCfilter}}, \code{\link{UCsmooth}}, 
#'          \code{\link{UCdisturb}}, \code{\link{UCcomponents}},
#'          \code{\link{UChp}}
#'          
#' @examples
#' y <- log(AirPassengers)
#' m1 <- UC(y)
#' m1 <- UC(y, model = "llt/different/arma(0,0)")
#' @rdname UC
#' @export
UC = function(y, u = NULL, model = "?/none/?/?", h = NA, outlier = NA, tTest = FALSE, criterion = "aic",
                   periods = NA, verbose = FALSE, stepwise = FALSE, p0 = NA, cLlik = TRUE, arma = TRUE){
    m1 = UCsetup(y, u, model, h, outlier, tTest, criterion, 
                 periods, verbose, stepwise, p0, cLlik, arma)
    m1 = UCestim(m1)
    m1 = UCvalidate(m1, verbose)
    m1 = UCdisturb(m1)
    m1 = UCsmooth(m1)
    m1 = UCcomponents(m1)
    return(m1)
}
