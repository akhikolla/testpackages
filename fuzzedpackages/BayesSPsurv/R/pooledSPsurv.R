#' @title pooledSPsurv
#' @description Markov Chain Monte Carlo (MCMC) to run Bayesian split population survival model with no frailties.
#'
#' @param duration survival stage equation written in a formula of the form Y ~ X1 + X2 + ... where Y is duration until failure or censoring.
#' @param immune split stage equation written in a formula of the form C ~ Z1 + Z2 + ... where C is a binary indicator of immunity.
#' @param Y0 the elapsed time since inception until the beginning of time period (t-1).
#' @param LY last observation year (coded as 1; 0 otherwise) due to censoring or failure.
#' @param data data.frame.
#' @param N number of MCMC iterations.
#' @param burn burn-in to be discarded.
#' @param thin thinning to prevent from autocorrelation.
#' @param w size of the slice in the slice sampling for (betas, gammas, rho). Write it as a vector. E.g. c(1,1,1).
#' @param m limit on steps in the slice sampling. A vector of values for beta, gamma, rho.
#' @param form type of parametric model (Weibull, Exponential, or Log-Logistic).
#'
#' @return pooledSPsurv returns an object of class \code{"SPsurv"}.
#'
#' A \code{"pooledSPsurv"} object has the following elements:
#' \item{betas}{matrix, numeric values of the posterior for each variable in the duration equation .}
#' \item{gammas}{matrix, numeric values of the posterior for each variable in the immune equation.}
#' \item{rho}{vector, numeric values of rho.}
#' \item{delta}{vector, numeric values of delta.}
#' \item{X}{matrix of X's variables.}
#' \item{Z}{matrix of Z's variables.}
#' \item{Y}{vector of `Y'.}
#' \item{Y0}{vector of `Y0'.}
#' \item{C}{vector of `C'.}
#' \item{form}{character, type of distribution.}
#' \item{call}{description for the model to be estimated.}
#'
#' @examples
#' \donttest{
#' walter <- spduration::add_duration(Walter_2015_JCR,"renewed_war",
#'                                    unitID = "ccode", tID = "year",
#'                                    freq = "year", ongoing = FALSE)
#'
#' set.seed(123456)
#'
#' model <-
#'     pooledSPsurv(
#'         duration = duration ~ fhcompor1 + lgdpl + comprehensive + victory +
#'             instabl + intensityln + ethfrac + unpko,
#'         immune   = cured ~ fhcompor1 + lgdpl + victory,
#'         Y0       = 't.0',
#'         LY       = 'lastyear',
#'         data     = walter,
#'         N        = 100,
#'         burn     = 10,
#'         thin     = 10,
#'         w        = c(1,1,1),
#'         m        = 10,
#'         form     = "Weibull"
#'     )
#'
#'
#' print(model)
#'
#' summary(model, parameter = "betas")
#'
#' # plot(model)
#'
#'
#'}
#'
#' @export

pooledSPsurv <- function(duration,
                  immune,
                  Y0,
                  LY,
                  data,
                  N,
                  burn,
                  thin,
                  w = c(1, 1, 1),
                  m = 10,
                  form = c('Weibull', 'exponential', 'loglog'))
{

    dis <- match.arg(form)
    model <- 'SPsurv'
    r   <- formcall(duration = duration, immune = immune, data = data, Y0 = Y0,
                    LY = LY, N = N, burn = burn, thin = thin, w = w,
                    m = m, form = dis, model = model)

    if(form == 'loglog') {
        results <- mcmcSPlog(Y = r$Y, Y0 = r$Y0, C = r$C, LY = r$LY, X = r$X, Z = r$Z,
                          N = r$N, burn = r$burn, thin = r$thin, w  = r$w, m  = r$m,
                          form = r$form)
    } else {
        results <- mcmcSP(Y = r$Y, Y0 = r$Y0, C = r$C, LY = r$LY, X = r$X, Z = r$Z,
                      N = r$N, burn = r$burn, thin = r$thin, w  = r$w, m  = r$m,
                      form = r$form)
    }

    results$call   <- match.call()
    class(results) <- c(class(results), model)
    results

}


#' @title summary.SPsurv
#' @description Returns a summary of a SPsurv object via \code{\link[coda]{summary.mcmc}}.
#' @param object an object of class \code{SPsurv}, the output of \code{\link{pooledSPsurv}}.
#' @param parameter one of three parameters of the \code{\link{pooledSPsurv}} output. Indicate either "betas," "gammas," or "lambda."
#' @param ... additional parameter
#' @return list. Empirical mean, standard deviation and quantiles for each variable.
#' @rdname pooledSPsurv
#' @export
#'
#'

summary.SPsurv <- function(object, parameter = c("betas", "gammas", "lambda"), ...){

    if (parameter == "betas")  sum <- summary(mcmc(object$betas),  ...)
    if (parameter == "gammas") sum <- summary(mcmc(object$gammas), ...)
    if (parameter == "lambda") sum <- summary(mcmc(object$lambda), ...)
    sum
}



#' @title print.SPsurv
#' @description Print method for a \code{\link{pooledSPsurv}} x.
#' @param x an x of class \code{SPsurv} (output of \code{\link{pooledSPsurv}}).
#' @rdname pooledSPsurv
#' @export

print.SPsurv <- function(x, ...){

    cat('Call:\n')
    print(x$call)
    cat('\n')
    x2 <- summary(x, parameter = 'betas')
    cat("\n", "Iterations = ", x2$start, ":", x2$end, "\n", sep = "")     # coda::mcmc
    cat("Thinning interval =", x2$thin, "\n")                             # coda::mcmc
    cat("Number of chains =", x2$nchain, "\n")                            # coda::mcmc
    cat("Sample size per chain =", (x2$end - x2$start)/x2$thin + 1, "\n") # coda::mcmc
    cat("\nEmpirical mean and standard deviation for each variable,")     # coda::mcmc
    cat("\nplus standard error of the mean:\n\n")
    cat('\n')
    cat('Duration equation: \n')
    print(summary(x, parameter = 'betas')$statistics)
    cat('\n')
    cat('Inmune equation: \n')
    print(summary(x, parameter = 'gammas')$statistics)
    cat('\n')

}



#' @title plot.SPsurv
#' @description Returns a plot of a pooledSPsurv object via \code{\link[coda]{plot.mcmc}}.
#' @param x an object of class \code{SPsurv}, the output of \code{\link{pooledSPsurv}}.
#' @param ... additional parameter.
#' @return list. Empirical mean, standard deviation and quantiles for each variable.
#' @rdname pooledSPsurv
#' @export
#'

plot.SPsurv <- function(x, ...){
    plot((coda::mcmc(x$betas)), ...)
}






