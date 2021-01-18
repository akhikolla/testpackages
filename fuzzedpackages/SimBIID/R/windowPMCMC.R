#' @title Time windows for \code{PMCMC} objects
#' 
#' @description \code{window} method for class \code{PMCMC}.
#' 
#' @details Acts as a wrapper function for \code{\link[coda]{window.mcmc}} 
#'          from the \code{coda} package
#' 
#' @param x a \code{PMCMC} object, usually as a result of a call to
#'          \code{PMCMC}.
#' @param \dots arguments to pass to \code{\link{window.mcmc}}
#' @return a \code{PMCMC} object
#'
#' @export
#' 
#' @seealso \code{\link{PMCMC}}, \code{\link{print.PMCMC}}, \code{\link{predict.PMCMC}}, \code{\link{summary.PMCMC}}
#'     \code{\link{plot.PMCMC}}
#' 
#' @method window PMCMC
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

window.PMCMC <- function(x, ...) {
    if(class(x) != "PMCMC"){
        stop("'x' is not a PMCMC object")
    }
    
    ## extract 'mcmc' object
    y <- x$pars
    
    ## extract subset
    y <- window(y, ...)
    
    ## generate new PMCMC object
    x$pars <- y
    x
}
