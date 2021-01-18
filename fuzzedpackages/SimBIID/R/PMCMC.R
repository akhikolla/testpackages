#' @title Runs particle MCMC algorithm
#'
#' @description Runs particle Markov chain Monte Carlo (PMCMC) algorithm using 
#'      a bootstrap particle filter for fitting infectious disease models to 
#'      time series count data.
#'
#' @details             Function runs a particle MCMC algorithm using a bootstrap particle filter for a given model. 
#'                      If running with \code{fixpars = TRUE} then this runs \code{niter} simulations
#'                      using fixed parameter values. This can be used to optimise the number of 
#'                      particles after a training run. Also has \code{print()}, \code{summary()},
#'                      \code{plot()}, \code{predict()} and \code{window()} methods.
#'
#' @param x 		    A \code{PMCMC} object, or a \code{data.frame} containing time series count data, 
#'                      with the first column called \code{t}, followed by columns of time-series counts. 
#'                      The time-series counts columns must be in the order of the `counts` object in the
#'                      `func` function (see below).
#' @param priors        A \code{data.frame} containing columns \code{parnames}, \code{dist}, \code{p1} and 
#'                      \code{p2}, with number of rows equal to the number of parameters. The column
#'                      \code{parname} simply gives names to each parameter for plotting and summarising.
#'                      Each entry in the \code{dist} column must contain one of \code{c("unif", "norm", "gamma")}, 
#'                      and the corresponding \code{p1} and \code{p2} entries relate to the hyperparameters 
#'                      (lower and upper bounds in the uniform case; mean and standard deviation in the 
#'                      normal case; and shape and rate in the gamma case).
#' @param func          A \code{SimBIID_model} object or an \code{XPtr} to simulation function. If the latter, 
#'                      then this function must take the following arguments in order: 
#'                      \itemize{
#'                      \item{\code{NumericVector pars}:}{ a vector of parameters;}
#'                      \item{\code{double tstart}:}{ the start time;}
#'                      \item{\code{double tstop}:}{ the end time;}
#'                      \item{\code{IntegerVector u}:}{ a vector of states at time \code{tstart};}
#'                      \item{\code{IntegerVector counts}:}{ a vector of observed counts at \code{tstop}.}}
#' @param u             A named vector of initial states.
#' @param npart         An integer specifying the number of particles for the bootstrap particle filter.
#' @param iniPars       A named vector of initial values for the parameters of the model. If left unspecified, 
#'                      then these are sampled from the prior distribution(s).
#' @param fixpars       A logical determining whether to fix the input parameters (useful for 
#'                      determining the variance of the marginal likelihood estimates).
#' @param niter         An integer specifying the number of iterations to run the MCMC.
#' @param nprintsum     Prints summary of MCMC to screen every \code{nprintsum} iterations. 
#'                      Also defines how often adaptive scaling of proposal variances occur.
#' @param adapt         Logical determining whether to use adaptive proposal or not.
#' @param propVar       A numeric (npars x npars) matrix with log (or logistic) covariances to use
#'                      as (initial) proposal matrix. If left unspecified then defaults to 
#'                      \code{diag(nrow(priors)) * (0.1 ^ 2) / nrow(priors)}.
#' @param adaptmixprop  Mixing proportion for adaptive proposal.
#' @param nupdate       Controls when to start adaptive update.
#' @param ...           Not used here.
#'
#' @return If the code throws an error, then it returns a missing value (\code{NA}). If 
#'         \code{fixpars = TRUE} it returns a list of length 2 containing:
#' \itemize{
#'      \item{\code{output}:}{ a matrix with two columns. The first contains the simulated
#'          log-likelihood, and the second is a binary indicator relating to whether the
#'          simulation was 'skipped' or not (1 = skipped, 0 = not skipped);}
#'      \item{\code{pars}:}{ a vector of parameters used for the simulations.}
#' }
#' If \code{fixpars = FALSE}, the routine returns a \code{PMCMC} object, essentially a 
#'          \code{list} containing:
#' \itemize{
#'  \item{\code{pars}:}{ an \code{mcmc} object containing posterior samples for the parameters;}
#'  \item{\code{u}:}{ a copy of the \code{u} input;}
#'  \item{\code{accrate}:}{ the cumulative acceptance rate;}
#'  \item{\code{npart}:}{ the chosen number of particles;}
#'  \item{\code{time}:}{ the time taken to run the routine (in seconds);}
#'  \item{\code{propVar}:}{ the proposal covariance for the parameter updates;}
#'  \item{\code{data}:}{ a copy of the \code{x} input;}
#'  \item{\code{priors}:}{ a copy of the \code{priors} input;}
#'  \item{\code{func}:}{ a copy of the \code{func} input.}
#' }
#' 
#' @rdname PMCMC
#' 
#' @references Andrieu C, Doucet A and Holenstein R (2010) <doi:10.1111/j.1467-9868.2009.00736.x>
#'
#' @export
#' 
#' @seealso \code{\link{print.PMCMC}}, \code{\link{plot.PMCMC}}, \code{\link{predict.PMCMC}}, \code{\link{summary.PMCMC}}
#'     \code{\link{window.PMCMC}}
#' 
#' @examples 
#' \donttest{
#' ## set up data to pass to PMCMC
#' flu_dat <- data.frame(
#'     t = 1:14,
#'     Robs = c(3, 8, 26, 76, 225, 298, 258, 233, 189, 128, 68, 29, 14, 4)
#' )
#' 
#' ## set up observation process
#' obs <- data.frame(
#'     dataNames = "Robs",
#'     dist = "pois",
#'     p1 = "R + 1e-5",
#'     p2 = NA,
#'     stringsAsFactors = FALSE
#' )
#' 
#' ## set up model (no need to specify tspan
#' ## argument as it is set in PMCMC())
#' transitions <- c(
#'     "S -> beta * S * I / (S + I + R + R1) -> I", 
#'     "I -> gamma * I -> R",
#'     "R -> gamma1 * R -> R1"
#' )
#' compartments <- c("S", "I", "R", "R1")
#' pars <- c("beta", "gamma", "gamma1")
#' model <- mparseRcpp(
#'     transitions = transitions, 
#'     compartments = compartments,
#'     pars = pars,
#'     obsProcess = obs
#' )
#' 
#' ## set priors
#' priors <- data.frame(
#'     parnames = c("beta", "gamma", "gamma1"), 
#'     dist = rep("unif", 3), 
#'     stringsAsFactors = FALSE)
#' priors$p1 <- c(0, 0, 0)
#' priors$p2 <- c(5, 5, 5)
#' 
#' ## define initial states
#' iniStates <- c(S = 762, I = 1, R = 0, R1 = 0)
#' 
#' set.seed(50)
#' 
#' ## run PMCMC algorithm
#' post <- PMCMC(
#'     x = flu_dat, 
#'     priors = priors,
#'     func = model, 
#'     u = iniStates,
#'     npart = 25,
#'     niter = 5000, 
#'     nprintsum = 1000
#' )
#' 
#' ## plot MCMC traces
#' plot(post, "trace")
#' 
#' ## continue for some more iterations
#' post <- PMCMC(post, niter = 5000, nprintsum = 1000)
#' 
#' ## plot traces and posteriors
#' plot(post, "trace")
#' plot(post)
#' 
#' ## remove burn-in
#' post <- window(post, start = 5000)
#' 
#' ## summarise posteriors
#' summary(post)
#' }
#' 

PMCMC <- function(x, ...) {
    UseMethod("PMCMC")
}

#' @rdname PMCMC
#' @export

PMCMC.PMCMC <- function(x, niter = 1000, nprintsum = 100, 
                        adapt = TRUE, adaptmixprop = 0.05, 
                        nupdate = 100, ...) {
    ## check object
    if(class(x) != "PMCMC") {
        stop("'x' not a PMCMC object")
    }
    if(length(x) <= 2){
        stop("'x' not a valid PMCMC object to continue running...")
    }
    
    ## collect arguments
    tempargs <- list(
        x = x$data, 
        priors = x$priors, 
        func = x$func, 
        u = x$u, 
        npart = x$npart,
        iniPars = x$pars[nrow(x$pars), -ncol(x$pars)], 
        fixpars = FALSE, 
        niter = niter, 
        nprintsum = nprintsum, 
        adapt = adapt, 
        propVar = x$propVar, 
        adaptmixprop = adaptmixprop, 
        nupdate = nupdate
    )
    
    ## run PMCMC
    temp <- do.call("PMCMC.default", tempargs)
    
    ## combine with original runs
    x$pars <- coda::as.mcmc(rbind(as.matrix(x$pars), as.matrix(temp$pars)))
    
    ## return new object
    x
}
    
#' @rdname PMCMC
#' @export

PMCMC.default <- function(
    x, priors, func, u, npart = 100,
    iniPars = NA, fixpars = FALSE, 
    niter = 1000, nprintsum = 100, 
    adapt = TRUE, propVar = NA, adaptmixprop = 0.05, nupdate = 100, ...
) {
    
    ## check inputs are present
    if(missing(x)){
        stop("'x' argument missing")
    }
    if(missing(priors)){
        stop("'priors' argument missing")
    }
    if(missing(func)){
        stop("'func' argument missing")
    }
    if(missing(u)){
        stop("'u' argument missing")
    } 
    
    ## check data set
    data <- x
    checkInput(data, "data.frame")
    if(colnames(data)[1] != "t"){
        stop("First column of 'data' must be 't'")
    }
    if(ncol(data) < 2) {
        stop("Must have at least one count column in 'data'")
    }
    checkInput(data$t, "numeric")
    for(j in 2:ncol(data)) {
        checkInput(data[, j, drop = TRUE], "numeric", int = TRUE)
    }  
    ## check time periods start at zero
    if(data$t[1] <= 0) {
        ## PROTECTS IN CASE INCIDENCE CURVES ARE SET IN SimBIID_model
        ## (since it resets each time that tstart > 0.0)
        stop("First time point in data set must be greater than zero for time being.")
    }
    ## append time zero to ensure bootstrap filter works
    data <- rbind(data[1, ], data)
    data[1, ] <- 0
    
    ## check priors 
    checkInput(priors, "data.frame", ncol = 4)
    if(!identical(colnames(priors), c("parnames", "dist", "p1", "p2"))) {
        stop("'priors' must have column names: 'parnames', 'dist', 'p1', 'p2'")
    }
    checkInput(priors$parnames, "character")
    checkInput(priors$dist, "character")
    checkInput(priors$p1, "numeric")
    checkInput(priors$p2, "numeric")
    checkInput(priors$dist, inSet = c("unif", "norm", "gamma"))
    temp <- priors[priors$dist == "unif", , drop = FALSE]
    if(nrow(temp) > 0) {
        ## check uniform bounds correct
        if(!all(apply(temp[, 3:4, drop = FALSE], 1, diff) > 0)) {
            stop("Priors: 'uniform' bounds not in correct order")
        }
    }
    temp <- priors[priors$dist == "norm", , drop = FALSE]
    if(nrow(temp) > 0) {
        ## check normal hyperparameters correct
        if(!all(temp$p2 > 0)){
            stop("Priors: 'normal' variance not positive")
        }
    }
    temp <- priors[priors$dist == "gamma", , drop = FALSE]
    if(nrow(temp) > 0) {
        ## check gamma bounds correct
        if(!all(temp$p1 > 0) | !all(temp$p2 > 0)){
            stop("Priors: 'gamma' hyperparameters not positive")
        }
    }
    
    ## check initial conditions
    if(!any(is.na(iniPars))) {
        checkInput(iniPars, c("vector", "numeric"), length = nrow(priors))
        if(!identical(names(iniPars), priors$parnames)){
            stop("'iniPars' names don't match to priors.")
        }
    } else {
        iniPars <- rep(NA, nrow(priors))
    }
    
    ## convert priors to matrix
    orig_priors <- priors
    priors$parnames <- NULL
    priors$dist <- match(priors$dist, c("unif", "norm", "gamma"))
    ## convert rate to scale for PMCMC
    priors$p2[priors$dist == "gamma"] <- 1 / priors$p2[priors$dist == "gamma"]
    priors <- as.matrix(priors)
    
    ## check function
    funcorig <- func
    if(class(func) != "XPtr" & class(func) != "SimBIID_model"){
        stop("'func' not a 'SimBIID_model' object or an 'XPtr' object")
    }
    
    ## check u
    checkInput(u, c("vector", "numeric"), int = TRUE, gte = 0)
    checkInput(sum(u), "numeric", int = TRUE, gt = 1)
    if(class(func) == "SimBIID_model") {
        checkInput(u, length = length(funcorig$compartments))
        if(!identical(names(u), funcorig$compartments)) {
            stop("'names(u)' does not match 'func$compartments'")
        }
    }
    
    if(class(func) == "XPtr"){
        RcppXPtrUtils::checkXPtr(func, "SEXP", c("NumericVector", "double", "double", "IntegerVector",
           "IntegerVector"))
    } else {
        ## check model
        if(func$tspan){
            warning("'SimBIID_model' object will have 'tspan' set to F\n")
        }
        if(is.null(func$obsProcess[1])){
            stop("'SimBIID_model' must have non-NULL 'obsProcess'")
        }
        if(!is.null(func$addVars[1])){
            stop("'SimBIID_model' can't have non-NULL 'addVars'")
        }
        if(!is.null(func$stopCrit[1])){
            stop("'SimBIID_model' can't have non-NULL 'stopCrit'")
        }
        if(!is.null(func$afterTstar[1])){
            stop("'SimBIID_model' can't have non-NULL 'afterTstar'")
        }
        ## check obsProcess against data
        if(!identical(func$obsProcess$dataNames, colnames(data)[-1])){
            stop("'data' columns do not match 'obsProcess$datanames'")
        }
        ## check priors against data
        if(!identical(orig_priors$parnames, func$pars)){
            stop("'pars' do not match 'priors'")
        }
        if(func$incidence) {
            func$compartments <- func$compartments[-grep("_inc$", func$compartments)]
        }
        ## generate model
        func <- mparseRcpp(
            transitions = func$transitions,
            compartments = func$compartments,
            pars = func$pars,
            obsProcess = func$obsProcess,
            addVars = NULL,
            stopCrit = NULL,
            tspan = FALSE,
            incidence = func$incidence,
            afterTstar = NULL,
            PF = TRUE,
            runFromR = FALSE
        )
        
        ## compile model
        func <- compileRcpp(func)
    }
    
    ## check proposal variances
    if(is.na(propVar[1])) {
        propVar <- diag(nrow(priors))
        ## adjust for scaling parameter for initial iterations
        propVar <- propVar * ((0.1 ^ 2) / nrow(propVar))
        propVar <- propVar / ((2.562 ^ 2) / nrow(propVar))
    } else {
        checkInput(propVar, c("numeric", "matrix"), nrow = nrow(priors), ncol = nrow(priors))
    }
    
    ## check runtime arguments
    checkInput(fixpars, c("logical", "vector"), 1)
    if(fixpars & any(is.na(iniPars))) {
        stop("Must input initial parameters if fixing parameters.")
    }
    checkInput(niter, c("numeric", "vector"), 1, int = TRUE, gt = 0)
    checkInput(npart, c("numeric", "vector"), 1, int = TRUE, gt = 0)
    checkInput(nprintsum, c("numeric", "vector"), 1, int = TRUE, gt = 0)
    
    ## check adaptive update and proposal covariance matrices
    checkInput(adapt, c("logical", "vector"), 1)
    checkInput(adaptmixprop, c("numeric", "vector"), 1, gt = 0, lt = 1)
    checkInput(nupdate, c("numeric", "vector"), 1, int = TRUE, gt = 0)
    
    ## run function
    output <- PMCMC_cpp(as.matrix(data), priors, orig_priors$parnames, iniPars, propVar, niter, npart, 
                    adaptmixprop, nprintsum, nupdate, as.numeric(fixpars), 
                    as.numeric(adapt), u, func);
    
    ## check to see if code has stopped afer initialisation
    if(length(output[[1]]) == 1) {
        return(NA)
    }
    if(length(output) == 2) {
        names(output) <- c("output", "pars")
        return(output)
    }
    
    ## convert output into correct format
    colnames(output[[1]]) <- c(orig_priors$parnames, "logPost")
    output[[1]] <- coda::as.mcmc(output[[1]])
    
    ## finalise output and set names
    output <- c(output[1], list(u), output[-1], list(data[-1, , drop = FALSE]), list(orig_priors), list(funcorig))
    names(output) <- c("pars", "u", "accrate", "npart", "time", "propVar", "data", "priors", "func")
        
    ## export class and object
    class(output) <- "PMCMC"
    output
}

