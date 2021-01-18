#' @title Summarises \code{ABCSMC} objects
#'
#' @description Summary method for \code{ABCSMC} objects.
#'
#' @param object    An \code{ABCSMC} object.
#' @param gen       The generation of ABC that you want to extract. If left missing then
#'                  defaults to final generation.
#' @param transfunc Is a \code{function} object where the arguments to the function must
#'                  match all or a subset of the parameters in the model. This function needs 
#'                  to return a \code{data.frame} object with columns containing the transformed
#'                  parameters.
#' @param ...           Not used here.
#'
#' @return          A \code{data.frame} with weighted posterior means and variances.
#'
#' @method summary ABCSMC
#'         
#' @export
#' 
#' @seealso \code{\link{ABCSMC}}, \code{\link{print.ABCSMC}}, \code{\link{plot.ABCSMC}}
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

summary.ABCSMC <- function(object, gen = NA, transfunc = NA, ...) {
    
    ## check x is an ABCSMC object
    if(class(object) != "ABCSMC"){
        stop("'object' not of type ABCSMC")
    }
    
    ## check gen is valid
    if(length(gen) != 1){
        stop("'gen' must be of length 1")
    }
    gen <- ifelse(is.na(gen), length(object$pars), gen)
    checkInput(gen, "numeric", 1, int = TRUE, gt = 0, lte = length(object$pars))
    
    ## extract relevant parts of the object
    weights <- object$weights[[gen]]
    pars <- object$pars[[gen]]
    pars <- as.data.frame(pars)
    
    ## check transformations
    if(length(transfunc) != 1){
        stop("'transfunc' must be of length 1")
    }
    if(is.function(transfunc)) {
    
        ## check function arguments
        fargs <- formals(transfunc)
        checkInput(names(fargs), inSet = colnames(pars))
        
        ## perform transformations if required
        temppars <- pars[, match(names(fargs), colnames(pars))]
        temppars <- as.data.frame(temppars)
        temppars <- as.list(temppars)
        names(temppars) <- names(fargs)
        temp <- do.call("transfunc", temppars)
        checkInput(temp, "data.frame", nrow = nrow(pars))
        
        ## bind to current posterior samples
        pars <- cbind(pars, temp)
    }
    
    ## extract parameter names
    parnames <- colnames(pars)
    
    ## calculate weighted mean
    postmn <- apply(cbind(weights, pars), 1, function(x) x[-1] * x[1])
    if(is.null(nrow(postmn))) postmn <- matrix(postmn, nrow = 1)
    postmn <- apply(postmn, 1, sum)
    
    ## calculate weighted variance
    postvar <- apply(rbind(postmn, pars), 2, function(x) (x[-1] - x[1])^2)
    postvar <- apply(cbind(weights, postvar), 1, function(x) x[-1] * x[1])
    if(is.null(nrow(postvar))) postvar <- matrix(postvar, nrow = 1)
    postvar <- apply(postvar, 1, sum)
    
    ## return summary
    postsum <- data.frame(Mean = postmn, SD = sqrt(postvar), ESS = object$ESS[[length(object$ESS)]])
    rownames(postsum) <- parnames
    postsum
}
    
    
    
