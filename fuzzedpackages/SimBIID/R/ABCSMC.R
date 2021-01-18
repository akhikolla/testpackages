#' @title Runs ABC-SMC algorithm
#'
#' @description Runs the Approximate Bayesian Computation Sequential Monte Carlo 
#'     (ABC-SMC) algorithm of Toni et al. (2009) for fitting 
#'     infectious disease models to time series count data. 
#'     
#' @details Samples initial particles from the specified prior distributions and 
#'     then runs a series of generations of ABC-SMC. The generations can either be 
#'     specified with a set of fixed tolerances, or by setting the tolerances at 
#'     each new generation as a quantile of the tolerances of the accepted particles 
#'     at the previous generation. Uses bisection method as detailed in McKinley et al. (2018).
#'     Passing an \code{ABCSMC} object into the \code{ABCSMC()} function acts as a 
#'     continuation method, allowing further generations to be run.
#'
#' @export
#'
#' @param x         An \code{ABCSMC} object or a named vector with entries
#'                  containing the observed summary statistics to match to. Names must match to `tols`.
#' @param priors    A \code{data.frame} containing columns \code{parnames}, \code{dist}, \code{p1} and 
#'                  \code{p2}, with number of rows equal to the number of parameters. The column
#'                  \code{parname} simply gives names to each parameter for plotting and summarising.
#'                  Each entry in the \code{dist} column must contain one of \code{c("unif", "norm", "gamma")}, 
#'                  and the corresponding \code{p1} and \code{p2} entries relate to the hyperparameters 
#'                  (lower and upper bounds in the uniform case; mean and standard deviation in the 
#'                  normal case; and shape and rate in the gamma case).
#' @param func      Function that runs the simulator and checks whether the simulation matches the data. 
#'                  The first four arguments must be \code{pars}, \code{data}, \code{tols} and 
#'                  \code{u}. If the simulations do not match the data then the function must 
#'                  return an \code{NA}, else it must returns a \code{vector} of simulated summary measures. 
#'                  In this latter case the output from the function must be a vector with length equal to 
#'                  \code{ncol(data)} and with entries in the same order as the columns of \code{data}.
#' @param u         A named vector of initial states.
#' @param tols 		A \code{vector} or \code{matrix} of tolerances, with the number of rows defining
#'                  the number of generations required, and columns defining the summary statistics
#'                  to match to. If a \code{vector}, then the length determines the summary statistics.
#'                  The columns/entries must match to those in `x`. 
#' @param ptols     The proportion of simulated outcomes at each generation to use to derive adaptive 
#'                  tolerances.
#' @param mintols   A vector of minimum tolerance levels.
#' @param ngen      The number of generations of ABC-SMC to run.
#' @param npart     An integer specifying the number of particles.
#' @param parallel  A \code{logical} determining whether to use parallel processing or not.
#' @param mc.cores  Number of cores to use if using parallel processing.
#' @param ...       Further arguments to pass to \code{func}. (Not used if extending runs.)
#'
#' @return An \code{ABCSMC} object, essentially a \code{list} containing:
#' \itemize{
#' \item{\code{pars}:}{ a \code{list} of \code{matrix} objects containing the accepted
#'                      particles. Each element of the list corresponds to a generation 
#'                      of ABC-SMC, with each matrix being of dimension 
#'                      \code{npart} x \code{npars};}
#' \item{\code{output}:}{ a \code{list} of \code{matrix} objects containing the simulated
#'                      summary statistics. Each element of the list corresponds to a
#'                      generation of ABC-SMC, with each matrix being of dimension 
#'                      \code{npart} x \code{ncol(data)};}
#' \item{\code{weights}:}{ a \code{list} of \code{vector} objects containing the particle
#'                      weights. Each element of the list corresponds to a
#'                      generation of ABC-SMC, with each vector being of length
#'                      \code{npart};}
#' \item{\code{ESS}:}{ a \code{list} of effective sample sizes. Each element of the list 
#'                      corresponds to a generation of ABC-SMC, with each vector being of 
#'                      length \code{npart};}
#' \item{\code{accrate}:}{ a \code{vector} of length \code{nrow(tols)} containing the
#'                      acceptance rates for each generation of ABC;}
#' \item{\code{tols}:}{ a copy of the \code{tols} input;}
#' \item{\code{ptols}:}{ a copy of the \code{ptols} input;}
#' \item{\code{mintols}:}{ a copy of the \code{mintols} input;}
#' \item{\code{priors}:}{ a copy of the \code{priors} input;}
#' \item{\code{data}:}{ a copy of the \code{data} input;}
#' \item{\code{func}:}{ a copy of the \code{func} input;}
#' \item{\code{u}}{ a copy of the \code{u} input;}
#' \item{\code{addargs}:}{ a copy of the \code{...} inputs.}
#' }
#' @rdname ABCSMC
#' 
#' @references Toni T, Welch D, Strelkowa N, Ipsen A and Stumpf MP (2009) <doi:10.1098/rsif.2008.0172>
#'     
#' @references McKinley TJ, Cook AR and Deardon R (2009) <doi:10.2202/1557-4679.1171>
#' 
#' @references McKinley TJ, Vernon I, Andrianakis I, McCreesh N, Oakley JE, Nsubuga RN, Goldstein M and White RG (2018) <doi:10.1214/17-STS618>
#'     
#' @seealso \code{\link{print.ABCSMC}}, \code{\link{plot.ABCSMC}}, \code{\link{summary.ABCSMC}}
#'     
#' @examples 
#' \donttest{
#' ## set up SIR simulationmodel
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
#' model <- compileRcpp(model)
#' 
#' ## generate function to run simulators
#' ## and return summary statistics
#' simSIR <- function(pars, data, tols, u, model) {
#' 
#'     ## run model
#'     sims <- model(pars, 0, data[2] + tols[2], u)
#'     
#'     ## this returns a vector of the form:
#'     ## completed (1/0), t, S, I, R (here)
#'     if(sims[1] == 0) {
#'         ## if simulation rejected
#'         return(NA)
#'     } else {
#'         ## extract finaltime and finalsize
#'         finaltime <- sims[2]
#'         finalsize <- sims[5]
#'     }
#'     
#'     ## return vector if match, else return NA
#'     if(all(abs(c(finalsize, finaltime) - data) <= tols)){
#'         return(c(finalsize, finaltime))
#'     } else {
#'         return(NA)
#'     }
#' }
#' 
#' ## set priors
#' priors <- data.frame(
#'     parnames = c("beta", "gamma"), 
#'     dist = rep("gamma", 2), 
#'     stringsAsFactors = FALSE
#' )
#' priors$p1 <- c(10, 10)
#' priors$p2 <- c(10^4, 10^2)
#' 
#' ## define the targeted summary statistics
#' data <- c(
#'     finalsize = 30, 
#'     finaltime = 76
#' )
#' 
#' ## set initial states (1 initial infection 
#' ## in population of 120)
#' iniStates <- c(S = 119, I = 1, R = 0)
#' 
#' ## set initial tolerances
#' tols <- c(
#'     finalsize = 50,
#'     finaltime = 50
#' )
#' 
#' ## run 2 generations of ABC-SMC
#' ## setting tolerance to be 50th
#' ## percentile of the accepted 
#' ## tolerances at each generation
#' post <- ABCSMC(
#'     x = data, 
#'     priors = priors, 
#'     func = simSIR, 
#'     u = iniStates, 
#'     tols = tols, 
#'     ptol = 0.2, 
#'     ngen = 2, 
#'     npart = 50,
#'     model = model
#' )
#' post
#' 
#' ## run one further generation
#' post <- ABCSMC(post, ptols = 0.5, ngen = 1)
#' post
#' summary(post)
#' 
#' ## plot posteriors
#' plot(post)
#' 
#' ## plot outputs
#' plot(post, "output")
#' }
#' 

ABCSMC <- function(x, ...) {
    UseMethod("ABCSMC")
}

#' @rdname ABCSMC
#' @export

ABCSMC.ABCSMC <- function(x, tols = NULL, ptols = NULL, mintols = NULL, ngen = 1, parallel = FALSE, mc.cores = NA, ...) {
    
    ## check inputs
    if(class(x) != "ABCSMC"){
        stop("'x' not ABCSMC object")
    }
    if(is.null(tols[1]) & is.null(ptols[1])){
        stop("Must input either 'tols' or 'ptols'")
    }
    if(!is.null(tols[1]) & !is.null(ptols[1])){
        stop("Must choose either 'tols' or 'ptols'")
    }
    
    ## set ne mintols if needed
    if(!is.null(mintols)) {
        x$mintols <- mintols
    }
    
    if(!is.null(tols[1])){
        ## extract tolerances and check against new tolerances
        if(!is.matrix(tols) & !is.vector(tols)){
            stop("'tols' not a vector or a matrix")
        }
        if(is.matrix(tols)) {
            checkInput(tols, ncol = length(x$data))
            if(!identical(colnames(tols), names(x$data))){
                stop("colnames(tols) does not match names(data)")
            }
        } else {
            checkInput(tols, length = length(x$data))
            if(!identical(names(tols), names(x$data))){
                stop("names(tols) does not match names(data)")
            }
        }
        checkInput(tols, c("numeric"))
        if(sum(apply(rbind(x$tols[nrow(x$tols), ], tols), 2, function(x) {
            sum(x[-1] > x[1])
        })) > 0) {
            stop("New tolerances not less than or equal to original tolerances")
        }
        ngen <- ifelse(is.null(nrow(tols)), 1, nrow(tols))
        temp <- rbind(x$tols[nrow(x$tols), ], tols)
        if(!all(apply(temp, 2, function(x) {
            all(diff(x) <= 0)
        }))){
            stop("'tols' cannot increase")
        }
        temp <- apply(temp, 2, function(x) {
            diff(x)
        })
        if(is.null(nrow(temp))) {
            temp <- matrix(temp, nrow = 1)
        }
        if(any(apply(temp, 1, function(x) !any(x < 0)))){
            stop("Some 'tols' must decrease between generations")
        }
    } else {
        checkInput(ngen, c("vector", "numeric"), 1, int = TRUE, gt = 0)
        checkInput(ptols, c("vector", "numeric"), 1, gt = 0, lt = 1)
        tols <- bisectTols(x$output[[length(x$output)]], x$data, x$tols[nrow(x$tols), ], ptols, ptollim = 0.1)
        names(tols) <- colnames(x$tols)
        tols <- ifelse(tols < 0, 0, tols)
        ## check against min tols
        if(!is.null(x$mintols[1])) {
            tols <- ifelse(tols < x$mintols, x$mintols, tols)
        }
        if(all(tols == x$tols[nrow(x$tols), ])){
            stop("Tolerances same as previous generation at this 'ptol'")
        }
    }
    
    ## collect arguments
    tempargs <- list(
        x = x$data, 
        npart = nrow(x$pars[[1]]), 
        tols = tols, 
        ptols = ptols,
        mintols = x$mintols,
        ngen = ngen,
        priors = x$priors, 
        func = x$func, 
        u = x$u,
        parallel = parallel,
        mc.cores = mc.cores, 
        prevPars = x$pars[[length(x$pars)]], 
        prevWeights = x$weights[[length(x$weights)]],
        prevESS = x$ESS[[length(x$ESS)]],
        genstart = nrow(x$tols) + 1
    )
    tempargs <- c(tempargs, x$addargs)
    
    ## run ABC-SMC
    temp <- do.call("ABCSMC.default", tempargs)
    
    ## combine with original runs
    x$pars <- c(x$pars, temp$pars)
    x$output <- c(x$output, temp$output)
    x$weights <- c(x$weights, temp$weights)
    x$ESS <- c(x$ESS, temp$ESS)
    x$accrate <- c(x$accrate, temp$accrate)
    x$tols <- rbind(x$tols, temp$tols)
    x$ptols <- temp$ptols
    
    ## return new object
    x
}

#' @rdname ABCSMC
#' @export

ABCSMC.default <- function(x, priors, func, u, tols = NULL, ptols = NULL, mintols = NULL,
                           ngen = 1, npart = 100, parallel = FALSE, mc.cores = NA, ...) {
    
    ## check missing arguments
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
    
    ## check inputs
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
    checkInput(npart, "numeric", 1, int = TRUE, gt = 1)
    checkInput(priors, "data.frame", ncol = 4)
    checkInput(func, "function", 1)
    data <- x
    checkInput(data, c("numeric", "vector"))
    if(is.null(names(data))){
        stop("'data' is not a named vector")
    }
    
    orig_tols <- tols
    if(is.null(tols[1]) & is.null(ptols[1])){
        stop("Must have at least one of 'tols' or 'ptols'")
    }
    if(!is.null(ptols[1])){
        checkInput(ngen, c("vector", "numeric"), 1, int = TRUE, gt = 0)
        checkInput(ptols, c("vector", "numeric"), 1, gt = 0, lt = 1)
        if(!is.vector(tols)){
            stop("'tols' must be a named vector if 'ptols' is set")
        }
    }
    if(!is.matrix(tols) & !is.vector(tols)){
       stop("'tols' not a matrix or vector")
    }
    if(is.matrix(tols)){
        checkInput(tols, "numeric", gte = 0, ncol = length(data))
        if(!identical(names(data), colnames(tols))){
            stop("colnames(tols) != names(data)")
        }
    } else {
        checkInput(tols, "numeric", gte = 0, length = length(data))
        if(!identical(names(data), names(tols))){
            stop("names(tols) != names(data)")
        }
        tols <- matrix(tols, nrow = 1)
        colnames(tols) <- names(data)
    }
    if(nrow(tols) > 1){
        if(!all(apply(tols, 2, function(x) {
            all(diff(x) <= 0)
        }))){
            stop("'tols' cannot increase")
        }
        temp <- apply(tols, 2, function(x) {
            diff(x)
        })
        if(is.null(nrow(temp))) {
            temp <- matrix(temp, nrow = 1)
        }
        if(apply(temp, 1, function(x) !any(x < 0))){
            stop("Some 'tols' must decrease between generations")
        }
    }
    if(is.null(ptols[1])){
        ngen <- nrow(tols)
    } else {
        if(nrow(tols) != 1){
            stop("Some weird error")
        }
        if(ngen > nrow(tols)) {
            tols <- rbind(tols, matrix(NA, ngen - 1, ncol(tols)))
        }
    }
    
    if(!is.null(mintols[1])) {
        checkInput(mintols, c("vector", "numeric"), gte = 0, length = ncol(tols))
        checkInput(names(mintols), c("vector", "character"), ident = colnames(tols))
        if(is.null(ptols[1])) {
            if(any(apply(rbind(mintols, tols), 2, function(x) {
                any(x[1] > x[-1])
            }))) {
                stop("'tols' can't be less than 'mintols'")
            }
        }
    }
    
    ## check function
    fargs <- formals(func)
    if(length(fargs) < 4){
        stop("Number of arguments of 'func' must be at least 4")
    }
    if(!identical(names(fargs)[1:4], c("pars", "data", "tols", "u"))){
        stop("First four arguments of 'func' must be: 'pars', 'data', 'tols' and 'u'")
    }
    ## check u
    checkInput(u, c("vector", "numeric"), int = TRUE, gte = 0)
    checkInput(sum(u), "numeric", int = TRUE, gt = 1)
    
    ## check priors
    if(!identical(colnames(priors), c("parnames", "dist", "p1", "p2"))){
        stop("Column names of 'priors' must be: 'parnames', 'dist', 'p1' and 'p2'")
    }
    checkInput(priors$parnames, "character")
    checkInput(priors$dist, "character")
    checkInput(priors$p1, "numeric")
    checkInput(priors$p2, "numeric")
    if(!all(priors$dist %in% c("unif", "norm", "gamma"))){
        stop("'priors' must be of form: 'unif', 'norm' or 'gamma'")
    }
    temp <- priors[priors$dist == "unif", , drop = FALSE]
    if(nrow(temp) > 0) {
        ## check uniform bounds correct
        if(!all(apply(temp[, 3:4, drop = FALSE], 1, diff) > 0)){
            stop("Priors: uniform bounds in wrong order")
        }
    }
    temp <- priors[priors$dist == "norm", , drop = FALSE]
    if(nrow(temp) > 0) {
        ## check normal hyperparameters correct
        if(!all(temp$p2 > 0)){
            stop("Priors: normal variances must be > 0")
        }
    }
    temp <- priors[priors$dist == "gamma", , drop = FALSE]
    if(nrow(temp) > 0) {
        ## check gamma bounds correct
        if(!all(temp$p1 > 0) | !all(temp$p2 > 0)){
            stop("Priors: gamma hyperparameters must be > 0")
        }
    }
    orig_priors <- priors
    priors$ddist <- paste0("d", priors$dist)
    priors$dist <- paste0("r", priors$dist)    
    
    ## set timer
    ptm_ini <- proc.time()

    ## set up output objects
    accrate <- rep(NA, ngen)
    pars <- list(); out <- list(); weights <- list();
    ESS <- list()
    
    ## extract ellipsis arguments
    args <- list(...)
    
    ## extract arguments for "func"
    fargs <- fargs[is.na(match(names(fargs), c("pars", "data", "tols", "u")))]
    if(length(fargs) > 0) {
        fargs1 <- match(names(fargs), names(args))
        if(!all(!is.na(fargs1))){
            stop(paste0("Need to include: '", paste(names(fargs)[is.na(fargs1)], collapse = "', '"), 
                        "' arguments in function call"))
        }
        fargs <- args[fargs1]
    }
    
    if(exists("prevPars", where = args)) {
        if(!exists("prevWeights", where = args)){
            stop("'prevWeights' does not exist")
        }
        if(!exists("prevESS", where = args)){
            stop("'prevESS' does not exist")
        }
        if(!exists("genstart", where = args)){
            stop("'genstart' does not exist")
        }
        ## set up dummy objects (removed at the end)
        pars[[1]] <- args$prevPars
        weights[[1]] <- args$prevWeights
        ESS[[1]] <- args$prevESS
        tols <- rbind(rep(NA, ncol(tols)), tols)
        out[[1]] <- NA
        propCov <- stats::cov(pars[[1]]) * 2
        init <- 2
        genstart <- args$genstart - 2
        ## remove additional arguments
        args <- args[-match("prevPars", names(args))]
        args <- args[-match("prevWeights", names(args))]
        args <- args[-match("prevESS", names(args))]
        args <- args[-match("genstart", names(args))]
    } else {
        init <- 1
        genstart <- 0
    }
    
    ## run sequential algorithm
    for(t in init:nrow(tols)) {  
        ## set up arguments
        if(t == 1) {
            tempWeights <- rep(1, npart)
            tempPars <- rep(1, npart)
        } else {
            tempWeights <- weights[[t - 1]]
            tempPars <- pars[[t - 1]]
        }
        runind <- T
        if(t != init){
            ## set tolerances if required
            if(!is.null(ptols[1])) {
                tols[t, ] <- bisectTols(out[[t - 1]], data, tols[t - 1, ], ptols, ptollim = 0.1)
                tols[t, ] <- ifelse(tols[t, ] < 0, 0, tols[t, ])
                ## check against min tols
                if(!is.null(mintols[1])) {
                    tols[t, ] <- ifelse(tols[t, ] < mintols, mintols, tols[t, ])
                }
                if(all(tols[t, ] == tols[t - 1, ])){
                    tols <- tols[1:(t - 1), , drop = FALSE]
                    runind <- F
                    message("Tolerances same as previous generation, so now stopping algorithm.")
                }
            }
        }
        
        ## set timer
        ptm <- proc.time()
        
        ## run generation
        if(runind){
            if(!parallel | npart == 1) {
                temp <- lapply(1:npart, runProp,
                    t = t, priors = priors, 
                    prevWeights = tempWeights, prevPars = tempPars, 
                    propCov = propCov, tols = tols[t, ], data = data, 
                    u = u,
                    func = func, func_args = fargs)
            } else  {
                ## set RNG generator to ensure reproducibility
                RNGkind("L'Ecuyer-CMRG")
                temp <- parallel::mclapply(1:npart, runProp,
                    t = t, priors = priors, 
                    prevWeights = tempWeights, prevPars = tempPars, 
                    propCov = propCov, tols = tols[t, ], data = data, 
                    u = u,
                    func = func, func_args = fargs, 
                    mc.cores = mc.cores)
                ## reset RNG generator
                RNGkind("default")
            }
            
            ## extract relative components
            weights[[t]] <- map_dbl(temp, "weightsNew")
            weights[[t]] <- weights[[t]] / sum(weights[[t]])
            ESS[[t]] <- 1 / sum(weights[[t]]^2)
            pars[[t]] <- map(temp, "pars")
            pars[[t]] <- do.call("rbind", pars[[t]])
            out[[t]] <- map(temp, "out")
            out[[t]] <- do.call("rbind", out[[t]])
            accrate[t] <- npart / sum(map_dbl(temp, "accrate"))
            
            ## set names
            colnames(pars[[t]]) <- priors$parnames
            colnames(out[[t]]) <- names(data)
            
            ## set proposal covariance
            propCov <- stats::cov(pars[[t]]) * 2
            
            ## stop timer
            ptm1 <- proc.time() - ptm
            
            ## print progress to the screen
            message(paste0("Generation ", t + genstart, ", accrate = ", signif(accrate[t], 2), 
                       ", time = ", signif(ptm1[3], 2), " secs\n"))
        }
    }
    
    ## stop timer
    ptm1 <- proc.time() - ptm_ini
    message(paste0("\nFinal run time = ", signif(ptm1[3], 2), " secs\n"))
    
    ## remove extraneous components if extending runs
    if(init > 1) {
        pars <- pars[-1]
        out <- out[-1]
        tols <- tols[-1, , drop = FALSE]
        weights <- weights[-1]
        ESS <- ESS[-1]
    }
    
    ## output results
    output <- list(pars = pars, output = out, weights = weights, ESS = ESS, accrate = accrate,
                   tols = tols, ptols = ptols, mintols = mintols, priors = orig_priors, data = data,
                   func = func, u = u, addargs = args)
    class(output) <- "ABCSMC"
    output
}

