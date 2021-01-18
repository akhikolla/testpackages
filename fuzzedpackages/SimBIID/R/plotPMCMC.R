#' @title Plots \code{PMCMC} objects
#'
#' @description Plot method for \code{PMCMC} objects.
#'
#' @param x             A \code{PMCMC} object.
#' @param type          Takes the value \code{"post"} if you want to plot posterior distributions.
#'                      Takes the value \code{"trace"} if you want to plot the trace plots.
#' @param joint         A logical describing whether joint or marginal distributions are wanted.
#' @param transfunc     Is a \code{function} object where the arguments to the function must
#'                      match all or a subset of the parameters in the model. This function needs 
#'                      to return a \code{data.frame} object with columns containing the transformed
#'                      parameters.
#' @param ask           Should the user ask before moving onto next trace plot.
#' @param ...           Not used here.
#'
#' @return A plot of the (approximate) posterior distributions obtained from fitting a particle 
#'      Markov chain Monte Carlo algorithm, or provides corresponding trace plots.
#'         
#' @method plot PMCMC
#' 
#' @export 
#' 
#' @seealso \code{\link{PMCMC}}, \code{\link{print.PMCMC}}, \code{\link{predict.PMCMC}}, \code{\link{summary.PMCMC}}
#'     \code{\link{window.PMCMC}}
#' 
#' @examples 
#' \donttest{
#' 
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

plot.PMCMC <- function(x, type = c("post", "trace"), joint = FALSE, transfunc = NA, ask = TRUE, ...) {
    
    ## check x
    if(class(x) != "PMCMC"){
        stop("'x' not a PMCMC object")
    }
    
    ## check type
    type <- type[1]
    checkInput(type, "character", inSet = c("post", "trace"))
    
    ## check joint
    checkInput(joint, c("vector", "logical"), 1)
    if(joint) {
        if(!requireNamespace("GGally", quietly = TRUE)) {
            stop("'GGally' package required for joint distribution plots")
        }
    }
    ## check ask
    checkInput(ask, c("vector", "logical"), 1)
    
    ## extract levels
    collev <- colnames(x$pars)
    
    if(type == "post") { 
        p <- as.matrix(x$pars) %>% as.data.frame()
                
        ## check for transformations if required
        checkInput(transfunc, length = 1, naAllow = TRUE)
        if(is.function(transfunc)) {
        
            ## check function arguments
            fargs <- formals(transfunc)
            checkInput(names(fargs), inSet = colnames(p))
            
            ## perform transformations if required
            temppars <- p[, match(names(fargs), colnames(p))]
            temppars <- as.list(temppars)
            names(temppars) <- names(fargs)
            temp <- do.call("transfunc", temppars)
            checkInput(temp, "data.frame", nrow = nrow(p))
            
            ## bind to current posterior samples
            p <- cbind(p, temp)
            collev <- colnames(p)
        } else {
            if(!is.na(transfunc)) {
                stop("'transfunc' incorrectly specified")
            }
        }
        
        ## plot posteriors
        if(!joint) {
             p <- p %>%
                gather(Parameter, value) %>%
                mutate(Parameter = factor(Parameter, levels = collev)) %>%
                ggplot(aes(x = value)) +
                    geom_density() +
                    xlab("Parameter value") + ylab("Density") +
                    facet_wrap(~ Parameter, scales = "free")
         } else {
            p <- GGally::ggpairs(p,
                diag = list(continuous = GGally::wrap("densityDiag", alpha = 0.8)),
                lower = list(continuous = "density"),
                upper = list(continuous = "blank"))
         }
         return(p)
     } else {
        plot(x$pars, density = FALSE, ask = ask)
     }
}
    
