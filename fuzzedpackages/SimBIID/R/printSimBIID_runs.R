#' @title Prints \code{SimBIID_runs} objects
#'
#' @description Print method for \code{SimBIID_runs} objects.
#'
#' @param x    A \code{SimBIID_runs} object.
#' @param ...           Not used here.
#'
#' @return Summary outputs printed to the screen.
#'  
#' @method print SimBIID_runs
#' 
#' @seealso \code{\link{mparseRcpp}}, \code{\link{plot.SimBIID_runs}}, \code{\link{run}}
#' 
#' @export
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
#'     pars = pars,
#'     tspan = TRUE
#' )
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
#' }
#' 

print.SimBIID_runs <- function(x, ...) {
    cat(paste("'SimBIID_runs' object with n =", nrow(x$sums), "replicates.\n"))
    if(nrow(x$sums) == 1){
        cat("\nOutput at final time point:\n")
        x$sums %>%
            dplyr::select(-rep) %>%
            as_tibble() %>%
            print()
        if(is.data.frame(x$runs)) {
            cat("\nTime-series counts:\n")
            x$runs %>%
                dplyr::select(-rep) %>%
                as_tibble() %>%
                print()
        }
    } else {
        cat("\nSummaries of outputs at final time point:\n")
        x$sums %>%
            dplyr::select(-rep) %>%
            summary() %>%
            print()
        # if(is.data.frame(x$runs)) {
        #     cat("\nSummaries of time-series counts:\n")
        #     x$runs %>%
        #         dplyr::select(-rep) %>%
        #         gather(output, value, -t) %>%
        #         group_by(t, output) %>%
        #         summarise(list(enframe(c(mean(value), quantile(value, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))))) %>%
        #         unnest() %>%
        #         spread(name, value) %>%
        #         print()
        # }
    }
}
