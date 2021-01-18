#' @title Predicts future course of outbreak from \code{PMCMC} objects
#'
#' @description Predict method for \code{PMCMC} objects.
#'
#' @param object             A \code{PMCMC} object.
#' @param tspan         A vector of times over which to output predictions.
#' @param npart         The number of particles to use in the bootstrap filter.
#' @param ...           Not used here.
#'
#' @return A \code{SimBIID_runs} object.
#'         
#' @method predict PMCMC
#' 
#' @export 
#' 
#' @seealso \code{\link{PMCMC}}, \code{\link{print.PMCMC}}, \code{\link{plot.PMCMC}}, \code{\link{summary.PMCMC}}
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
#' ## run PMCMC algorithm for first three days of data
#' post <- PMCMC(
#'     x = flu_dat[1:3, ], 
#'     priors = priors, 
#'     func = model, 
#'     u = iniStates, 
#'     npart = 75, 
#'     niter = 10000, 
#'     nprintsum = 1000
#' )
#' 
#' ## plot traces
#' plot(post, "trace")
#' 
#' ## run predictions forward in time
#' post_pred <- predict(
#'     window(post, start = 2000, thin = 8), 
#'     tspan = 4:14
#' )
#' 
#' ## plot predictions
#' plot(post_pred, quant = c(0.6, 0.75, 0.95))
#' }
#' 

predict.PMCMC <- function(object, tspan, npart = 50, ...) {
    
    ## check object
    if(class(object) != "PMCMC"){
        stop("'object' not a PMCMC object")
    }
    
    if(class(object$func) != "SimBIID_model"){
        stop("'object$func' not a SimBIID_model object")
    }
    
    if(missing(tspan)){
        stop("'tspan' can't be missing")
    }
    
    ## check stop
    checkInput(tspan, c("vector", "numeric"), gt = max(object$data$t))
    tspan <- sort(tspan)
    
    ## check npart
    checkInput(npart, c("vector", "numeric"), 1, int = TRUE, gt = 0)
    
    ## extract parameters and remove extraneous columns
    pars <- as.matrix(object$pars) %>%
        as.data.frame() %>%
        select(one_of(object$func$pars)) %>%
        as.matrix()
    
    ## run bootstrap filter to get states at each time point of the dataset
    prevStates <- bootStates(object$data, object$func, pars, object$u, npart)
    pars <- prevStates$pars
    prevStates <- prevStates$output
    
    ## extract final states
    iniStates <- prevStates %>%
        group_by(rep) %>%
        slice(n()) %>%
        ungroup() %>%
        select(one_of(object$func$compartments)) %>%
        as.matrix()
    
    ## generate model
    func <- object$func
    if(!func$tspan){
        message("For predictions 'SimBIID_model' object will have 'tspan' set to T\n")
    }
    # if(!is.null(func$obsProcess[1])){
    #     message("For predictions, 'obsProcess' will be removed from 'SimBIID_model' simulations\n")
    # }
    if(!is.null(func$addVars[1])){
        stop("'SimBIID_model' can't have non-NULL 'addVars'")
    }
    if(!is.null(func$stopCrit[1])){
        stop("'SimBIID_model' can't have non-NULL 'stopCrit'")
    }
    if(func$incidence) {
        func$compartments <- func$compartments[-grep("_inc$", func$compartments)]
    }
    func <- mparseRcpp(
        transitions = func$transitions,
        compartments = func$compartments,
        pars = func$pars,
        obsProcess = func$obsProcess,
        addVars = NULL,
        stopCrit = NULL,
        tspan = TRUE,
        incidence = func$incidence,
        afterTstar = NULL,
        PF = FALSE,
        runFromR = TRUE
    )
    compfunc <- compileRcpp(func)
    
    ## use run method to produce forward predictions
    out <- list()
    outsums <- list()
    for(i in 1:nrow(pars)) {
        out[[i]] <- compfunc(pars[i, ], max(object$data$t), max(tspan), iniStates[i, ], tspan)
        outsums[[i]] <- out[[i]][[1]]
        out[[i]] <- out[[i]][[2]]
        outsums[[i]] <- as.data.frame(matrix(outsums[[i]], nrow = 1))
        out[[i]] <- as.data.frame(out[[i]])
        tempnames <- c("completed", "t", colnames(iniStates))
        if(is.data.frame(func$obsProcess)) {
            tempnames <- c(tempnames, func$obsProcess$dataNames)
        }
        colnames(outsums[[i]]) <- tempnames
        tempnames <- c("t", colnames(iniStates))
        if(is.data.frame(func$obsProcess)) {
            tempnames <- c(tempnames, func$obsProcess$dataNames)
        }
        colnames(out[[i]]) <- tempnames
    }
    ## bind to prevStates
    out <- out %>%
        bind_rows(.id = "rep") %>%
        rbind(prevStates) %>%
        arrange(rep, t)
    
    ## extract summaries
    outsums <- outsums %>%
        bind_rows(.id = "rep")
    
    ## bind as output list and return
    out <- list(sums = outsums, runs = out, bootEnd = max(object$data$t))
    class(out) <- "SimBIID_runs"
    out
}
    
