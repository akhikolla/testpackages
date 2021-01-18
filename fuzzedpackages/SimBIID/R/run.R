#' @title Runs \code{SimBIID_model} object
#'
#' @description Wrapper function that compiles (if necessary) and runs
#'              a \code{SimBIID_model} object. Returns results in a 
#'              user-friendly manner as a \code{SimBIID_run} object, 
#'              for which \code{print()} and \code{plot()} generics 
#'              are provided.
#'
#' @param model An object of class \code{SimBIID_model}.
#' @param pars A named vector of parameters.
#' @param tstart The time at which to start the simulation.
#' @param tstop The time at which to stop the simulation.
#' @param u A named vector of initial states.
#' @param tspan A numeric vector containing the times at which to 
#'               save the states of the system.
#' @param nrep Specifies the number of simulations to run.
#' @param parallel  A \code{logical} determining whether to use parallel processing or not.
#' @param mc.cores  Number of cores to use if using parallel processing.
#'
#' @return An object of class \code{SimBIID_run}, essentially a list 
#'         containing elements:
#'         \itemize{
#'             \item{sums:}{ a \code{data.frame()} with summaries of the model runs. This
#'             includes columns \code{run}, \code{completed}, \code{t}, \code{u*} 
#'             (see help file for \code{SimBIID_model} for more details);}
#'             \item{runs:}{ a \code{data.frame()} object, containing columns: \code{run},
#'             \code{t}, \code{u*} (see help file for \code{SimBIID_model} for more details).
#'             These contain time series counts for the simulations. Note that this will
#'             only be returned if \code{tspan = TRUE} in the original \code{SimBIID_model} object.}
#'             \item{bootEnd:}{ a time point denoting when bootstrapped estimates end and predictions
#'             begin (for \code{predict.PMCMC()} method). }
#'         } 
#' 
#' @export
#' 
#' @seealso \code{\link{mparseRcpp}}, \code{\link{print.SimBIID_runs}}, \code{\link{plot.SimBIID_runs}}
#' 
#' @examples 
#' \donttest{
#' ## set up SIR simulation model
#' transitions <- c(
#'     "S -> beta * S * I -> I", 
#'     "I -> gamma * I -> R"
#' )
#' compartments <- c("S", "I", "R")
#' pars <- c("beta", "gamma")
#' model <- mparseRcpp(
#'     transitions = transitions, 
#'     compartments = compartments,
#'     pars = pars
#' )
#' 
#' ## compile and run model
#' sims <- run(
#'     model = model,
#'     pars = c(beta = 0.001, gamma = 0.1),
#'     tstart = 0,
#'     tstop = 100,
#'     u = c(S = 119, I = 1, R = 0)
#' )
#' sims
#' 
#' ## add tspan option to return
#' ## time series counts at different
#' ## time points
#' model <- mparseRcpp(
#'     transitions = transitions, 
#'     compartments = compartments,
#'     pars = pars,
#'     tspan = TRUE
#' )
#' sims <- run(
#'     model = model,
#'     pars = c(beta = 0.001, gamma = 0.1),
#'     tstart = 0,
#'     tstop = 100,
#'     u = c(S = 119, I = 1, R = 0),
#'     tspan = seq(1, 100, length.out = 10)
#' )
#' sims
#' 
#' ## run 100 replicate simulations and
#' ## plot outputs
#' sims <- run(
#'     model = model,
#'     pars = c(beta = 0.001, gamma = 0.1),
#'     tstart = 0,
#'     tstop = 100,
#'     u = c(S = 119, I = 1, R = 0),
#'     tspan = seq(1, 100, length.out = 10),
#'     nrep = 100
#' )
#' sims
#' plot(sims)
#' }

run <- function(
    model,
    pars,
    tstart,
    tstop,
    u,
    tspan,
    nrep = 1,
    parallel = FALSE,
    mc.cores = NA
) {
    ## check inputs
    if(missing(model)) {
        stop("'model' object missing")
    }
    if(class(model) != "SimBIID_model") {
        stop("'model' object not of class 'SimBIID_model'")
    }
    if(missing(pars)) {
        stop("'pars' object missing")
    }
    if(missing(tstart)) {
        stop("'tstart' object missing")
    }
    if(missing(tstop)) {
        stop("'tstop' object missing")
    }
    if(missing(u)) {
        stop("'u' object missing")
    }
    ## check inputs
    if(!model$runFromR) {
        stop("'model' must be specified with 'runFromR = TRUE'")
    }
    checkInput(pars, c("vector", "numeric"), length(model$pars))
    parnames <- names(pars)
    if(!identical(parnames, model$pars)){
        stop("Names of 'pars' do not match model specification")
    }
    checkInput(tstart, c("vector", "numeric"), 1)
    checkInput(tstop, c("vector", "numeric"), 1, gt = tstart)
    checkInput(u, c("vector", "numeric"), length(model$compartments), int = TRUE, gte = 0)
    if(sum(u) <= 0){
        stop("'sum(u)' not greater than zero.")
    }
    unames <- names(u)
    if(!identical(unames, model$compartments)){
        stop("Names of 'u' do not match model specification")
    }
    if(missing(tspan) & model$tspan){
        stop("'SimBIID_model' requires that 'tspan' must be set")
    }
    if(!missing(tspan) & !model$tspan){
        stop("'SimBIID_model' does not allow a 'tspan' argument")
    }
    if(!missing(tspan)){
        checkInput(tspan, c("vector", "numeric"), int = TRUE, gte = tstart)
        tspan <- sort(tspan)
    }
    checkInput(nrep, "numeric", 1, int = TRUE, gt = 0)
    checkInput(parallel, c("vector", "logical"), 1)
    if(parallel) {
        if(!requireNamespace("parallel", quietly = TRUE)) {
            stop("Must have 'parallel' package installed to use parallelisation")
        }
        nc <- parallel::detectCores()
        nc <- ifelse(is.na(nc), 1, nc)
        if(!is.na(mc.cores[1])) {
            checkInput(mc.cores, "numeric", 1, int = TRUE)
            mc.cores <- min(nc, mc.cores)
        } else {
            mc.cores <- nc
        }
        parallel <- (mc.cores > 1)
        message(paste0("Number of cores: ", mc.cores, "\n"))
    }
    
    ## compile model
    compModel <- compileRcpp(model)
    
    ## run simulations
    if(!parallel | nrep == 1) {
        if(!missing(tspan)) {
            sims <- lapply(1:nrep, function(i, model, pars, tstart, tstop, u, tspan){
                model(pars, tstart, tstop, u, tspan)
            }, model = compModel, pars = pars, tstart = tstart, tstop = tstop, u = u, tspan = tspan)
        } else {
            sims <- lapply(1:nrep, function(i, model, pars, tstart, tstop, u){
                model(pars, tstart, tstop, u)
            }, model = compModel, pars = pars, tstart = tstart, tstop = tstop, u = u)
        }
    } else {
        ## set RNG generator to ensure reproducibility
        RNGkind("L'Ecuyer-CMRG")
        if(!missing(tspan)) {
            sims <- parallel::mclapply(1:nrep, function(i, model, pars, tstart, tstop, u, tspan){
                model(pars, tstart, tstop, u, tspan)
            }, model = compModel, pars = pars, tstart = tstart, tstop = tstop, u = u, tspan = tspan, mc.cores = mc.cores)
        } else {
            sims <- parallel::mclapply(1:nrep, function(i, model, pars, tstart, tstop, u){
                model(pars, tstart, tstop, u)
            }, model = compModel, pars = pars, tstart = tstart, tstop = tstop, u = u, mc.cores = mc.cores)
        }
        ## reset RNG generator
        RNGkind("default")
    }
    ## extract simulations in useful format
    if(missing(tspan)) {
        sims <- do.call("rbind", sims)
        sims <- cbind(1:nrow(sims), sims)
        tempnames <- c("rep", "completed", "t", model$compartments)
        if(is.data.frame(model$obsProcess)) {
            tempnames <- c(tempnames, model$obsProcess$dataNames)
        }
        colnames(sims) <- tempnames
        sums <- as.data.frame(sims)
        runs <- NA
    } else {
        sums <- do.call("rbind", lapply(sims, function(x){
            x[[1]]
        }))
        sums <- cbind(1:nrow(sums), sums)
        tempnames <- c("rep", "completed", "t", model$compartments)
        if(is.data.frame(model$obsProcess)) {
            tempnames <- c(tempnames, model$obsProcess$dataNames)
        }
        colnames(sums) <- tempnames
        sums <- as.data.frame(sums)
        
        sims <- do.call("rbind", lapply(1:length(sims), function(i, x){
            x <- x[[i]][[2]]
            x <- cbind(rep(i, nrow(x)), x)
            x
        }, x = sims))
        tempnames <- c("rep", "t", model$compartments)
        if(is.data.frame(model$obsProcess)) {
            tempnames <- c(tempnames, model$obsProcess$dataNames)
        }
        colnames(sims) <- tempnames
        runs <- as.data.frame(sims)
    }
    ## return list of simulated outputs
    out <- list(sums = sums, runs = runs, bootEnd = NA)
    class(out) <- "SimBIID_runs"
    out
}
