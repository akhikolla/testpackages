#' @title Produces ABC reference table
#'
#' @description Produces reference table of simulated outcomes for use in various
#'     Approximate Bayesian Computation (ABC) algorithms.
#' 
#' @details Runs simulations for a large number of particles, either pre-specified or
#'     sampled from the a set of given prior distributions. Returns a table of summary 
#'     statistics for each particle. Useful for deciding on initial tolerances during an
#'     \code{\link{ABCSMC}} run, or for producing a reference table to use in e.g. the 
#'     ABC with Random Forests approach of Raynal et al. (2017).
#'
#' @export
#'
#' @param npart     The number of particles (must be a positive integer).
#' @param priors    A \code{data.frame} containing columns \code{parnames}, \code{dist}, \code{p1} and 
#'                  \code{p2}, with number of rows equal to the number of parameters. The column
#'                  \code{parname} simply gives names to each parameter for plotting and summarising.
#'                  Each entry in the \code{dist} column must contain one of \code{c("unif", "norm", "gamma")}, 
#'                  and the corresponding \code{p1} and \code{p2} entries relate to the hyperparameters 
#'                  (lower and upper bounds in the uniform case; mean and standard deviation in the 
#'                  normal case; and shape and rate in the gamma case).
#' @param pars      A named vector or matrix of parameters to use for the simulations. If \code{pars} is a vector then
#'                  this is repeated `npart` times, else it must be a \code{matrix} with `npart` rows. 
#'                  You cannot specify both `pars` and `priors`.
#' @param func      Function that runs the simulator. The first argument must be \code{pars}. The function
#'                  must return a \code{vector} of simulated summary measures, or a missing value (\code{NA})
#'                  if there is an error. The output from the function must be a vector with length equal 
#'                  to \code{length(sumNames)}.
#' @param sumNames  A \code{character} vector of summary statistic names.
#' @param parallel  A \code{logical} determining whether to use parallel processing or not.
#' @param mc.cores  Number of cores to use if using parallel processing.
#' @param ...       Extra arguments to be passed to \code{func}.
#'
#' @return An \code{data.frame} object with \code{npart} rows, where the first \code{p} columns correspond to 
#'         the proposed parameters, and the remaining columns correspond to the simulated outputs.
#'         
#' @references Raynal, L, Marin J-M, Pudlo P, Ribatet M, Robert CP and Estoup A. (2017) <ArXiv:1605.05537>
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
#' model <- compileRcpp(model)
#' 
#' ## generate function to run simulators
#' ## and produce final epidemic size and time
#' ## summary statistics
#' simRef <- function(pars, model) {
#'     ## run model over a 100 day period with
#'     ## one initial infective in a population
#'     ## of 120 individuals
#'     sims <- model(pars, 0, 100, c(119, 1, 0))
#'     
#'     ## return vector of summary statistics
#'     c(finaltime = sims[2], finalsize = sims[5])
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
#' ## produce reference table by sampling from priors
#' ## (add additional arguments to 'func' at the end)
#' refTable <- ABCRef(
#'     npart = 100, 
#'     priors = priors, 
#'     func = simRef, 
#'     sumNames = c("finaltime", "finalsize"),
#'     model = model
#' )
#' refTable
#' }
#' 

ABCRef <- function(npart, priors, pars, func, sumNames, parallel = FALSE, mc.cores = NA, ...) {
    
    ## check inputs
    if(missing(npart)){
        stop("'npart' must be provided")
    }
    if(missing(priors) & missing(pars)){
        stop("'priors' or 'pars' must be provided")
    }
    if(missing(func)){
        stop("'func' must be provided")
    }
    if(missing(sumNames)){
        stop("'sumNames' must be provided")
    }
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
    checkInput(npart, "numeric", 1, int = TRUE)
    if(!missing(priors)){
        if(!missing(pars)) {
            stop("Don't provide 'priors' and 'pars', you must choose one or the other.")
        }
        checkInput(priors, "data.frame", ncol = 4)
    } else {
        if(!is.vector(pars) & !is.matrix(pars)){
            stop("'pars' is neither a vector or a matrix")
        }
        checkInput(pars, "numeric")
        if(is.matrix(pars)){
            if(nrow(pars) != 1 & nrow(pars) != npart) {
                stop("'pars' must have either 1 row or 'npart' rows, or be a vector")
            }
            if(is.null(colnames(pars))){
                stop("'pars' must have column names")
            }
        } else {
            if(is.null(names(pars))){
                stop("'pars' must be named")
            }
        }
    }
    checkInput(func, "function", 1)
    checkInput(sumNames, c("vector", "character"))
    fargs <- formals(func)
    if(length(fargs) < 1) {
        stop("'func' must contain more than one argument")
    }
    if(names(fargs)[1] != "pars"){
        stop("First argument to 'func' must be 'pars'")
    }
    checkInput(npart, gt = 1)
    
    ## check priors
    if(!missing(priors)){
        if(!identical(colnames(priors), c("parnames", "dist", "p1", "p2"))){
            stop("colnames(priors) must be: 'parnames', 'dist', 'p1' and 'p2'")
        }
        checkInput(priors$parnames, "character")
        checkInput(priors$dist, "character")
        checkInput(priors$p1, "numeric")
        checkInput(priors$p2, "numeric")
        checkInput(priors$dist, inSet = c("unif", "norm", "gamma"))
        temp <- priors[priors$dist == "unif", , drop = FALSE]
        if(nrow(temp) > 0) {
            ## check uniform bounds correct
            if(!all(apply(temp[, 3:4, drop = FALSE], 1, diff) > 0)){
                stop("Priors: uniform bounds in incorrect order")
            }
        }
        temp <- priors[priors$dist == "norm", , drop = FALSE]
        if(nrow(temp) > 0) {
            ## check normal hyperparameters correct
            if(!all(temp$p2 > 0)){
                stop("Priors: normal variance must be > 0")
            }
        }
        temp <- priors[priors$dist == "gamma", , drop = FALSE]
        if(nrow(temp) > 0) {
            ## check gamma bounds correct
            if(!all(temp$p1 > 0) | !all(temp$p2 > 0)){
                stop("Priors: gamma hyperparameters must be > 0")
            }
        }
        priors$ddist <- paste0("d", priors$dist)
        priors$dist <- paste0("r", priors$dist)
        parnames <- priors$parnames
    } else {
        parnames <- colnames(pars)
        if(is.matrix(pars)){
            if(nrow(pars) == 1){
                pars <- as.vector(pars)
            }
        }
        if(is.vector(pars)){
            priors <- matrix(rep(pars, npart), ncol = length(pars), byrow = TRUE)
            priors <- as.data.frame(priors)
        } else {
            priors <- pars
        }
    }
    
    ## extract arguments for "func"
    fargs <- fargs[is.na(match(names(fargs), "pars"))]
    if(length(fargs) > 0) {
        args <- list(...)
        fargs1 <- match(names(fargs), names(args))
        if(!all(!is.na(fargs1))){
            stop(paste0("Need to include: ", paste(names(fargs)[is.na(fargs1)], collapse = ", "), 
                        " arguments in function call"))
        }
        fargs <- args[fargs1]
    }
    
    ## set timer
    ptm_ini <- proc.time()
    ptm <- ptm_ini

    ## set up output objects
    pars <- list(); out <- list(); 
    
    ## run generation
    if(!parallel | npart == 1) {
        temp <- lapply(1:npart, runRef, priors = priors, func = func, func_args = fargs)
    } else  {
        ## set RNG generator to ensure reproducibility
        RNGkind("L'Ecuyer-CMRG")
        temp <- parallel::mclapply(1:npart, runRef, priors = priors, func = func, func_args = fargs, 
                         mc.cores = mc.cores)
        ## reset RNG generator
        RNGkind("default")
    }
    
    ## extract relative components
    out <- map(temp, "out")
    pars <- map(temp, "pars")
    pars <- pars[!map_lgl(out, function(x) is.na(x[1]))]
    out <- out[!map_lgl(out, function(x) is.na(x[1]))]
    
    while(length(out) != npart) {
        
        ## stop timer
        ptm1 <- proc.time() - ptm
        ptm <- ptm1
        
        ## print progress
        message(paste0("Current pars = ", length(out), ", time = ", signif(ptm1[3], 2), " secs\n"))
        
        ## run generation
        if(!parallel) {
            temp <- lapply(1:npart, runRef, priors = priors, func = func)
        } else  {
            temp <- parallel::mclapply(1:npart, runRef, priors = priors, func = func, mc.cores = mc.cores)
        }
        ## extract relative components
        out <- c(out, map(temp, "out"))
        pars <- c(pars, map(temp, "pars"))
        pars <- pars[!map_lgl(out, function(x) is.na(x[1]))]
        out <- out[!map_lgl(out, function(x) is.na(x[1]))]
    }
    out <- do.call("rbind", out)
    pars <- do.call("rbind", pars)
    out <- as.data.frame(out)
    pars <- as.data.frame(pars)
    
    ## set names
    colnames(pars) <- parnames
    if(ncol(out) != length(sumNames)){
        stop("ncol(out) != length(sumNames)")
    }
    colnames(out) <- sumNames
    
    ## stop timer
    ptm1 <- proc.time() - ptm_ini
    
    ## print progress to the screen
    message(paste0("Final run time = ", signif(ptm1[3], 2), " secs\n"))
    
    return(cbind(pars, out))
}

