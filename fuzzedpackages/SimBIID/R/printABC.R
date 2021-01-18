#' @title Prints \code{ABCSMC} objects
#'
#' @description Print method for \code{ABCSMC} objects.
#'
#' @param x    An \code{ABCSMC} object.
#' @param ...           Not used here.
#'
#' @return Summary outputs printed to the screen.
#' 
#' @method print ABCSMC
#'         
#' @export
#' 
#' @seealso \code{\link{ABCSMC}}, \code{\link{plot.ABCSMC}}, \code{\link{summary.ABCSMC}}
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

print.ABCSMC <- function(x, ...) {
    ## check object is a ABCSMC object
    if(class(x) != "ABCSMC"){
        stop("'x' is not a ABCSMC object")
    }
    
    cat("An object of class: 'ABCSMC'\n")
    cat(paste0("Consists of ", nrow(x$tols), " generations with ", nrow(x$priors), " parameters.\n"))
    
    ## print data information
    cat("\nData:\n")
    print(x$data, row.names = FALSE)
    
    ## print tolerance information
    temp <- x$tols %>%
        as.data.frame() %>%
        mutate(Generation = 1:n()) %>%
        dplyr::select(Generation, everything()) %>%
        mutate(ESS = do.call("c", x$ESS))
    cat("\nTolerances:\n")
    print(temp, row.names = FALSE)
    
    ## print prior information
    temp <- x$priors %>%
        mutate(p1 = as.character(signif(p1, 2))) %>%
        mutate(p2 = as.character(signif(p2, 2))) %>%
        mutate(temp = ifelse(dist == "unif", paste0("U(lower = ", p1, ", upper = ", p2, ")"), NA)) %>%
        mutate(temp = ifelse(dist == "gamma", paste0("G(shape = ", p1, ", rate = ", p2, ")"), temp)) %>%
        mutate(temp = ifelse(dist == "norm", paste0("N(mean = ", p1, ", sd = ", p2, ")"), temp)) %>%
        mutate(temp = paste0(parnames, " ~ ", temp)) %>%
        dplyr::select(temp)
    colnames(temp) <- ""
    cat("\nPriors:\n")
    print(temp, row.names = FALSE, col.names = FALSE, quote = FALSE)
}

