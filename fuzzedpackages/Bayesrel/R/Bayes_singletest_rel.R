#'
#' calculate single test reliability estimates
#' @description calculate Bayesian and frequentist single test reliability measures.
#' Reported are Bayesian credible intervals (HDI) and frequentist confidence intervals (non parametric or parametric bootstrap).
#' The estimates supported are Cronbach alpha, lambda2/4/6, the glb, and Mcdonald omega. Beware of lambda4 with many indicators,
#' the computational effort is considerable
#'
#' @param data The dataset to be analyzed, observations are rows, items are columns
#' @param estimates A character vector containing the estimands, we recommend using lambda4 with only a few items
#' due to the computation time
#' @param cov.mat A covariance matrix can be supplied instead of a dataset,
#' but number of observations needs to be specified
#' @param interval A number specifying the uncertainty interval
#' @param n.iter A number for the iterations of the Gibbs Sampler
#' @param n.burnin A number for the burnin in the Gibbs Sampler
#' @param thin A number for the thinning of the MCMC samples
#' @param n.chains A number for the chains to run for the MCMC sampling
#' @param n.boot A number for the bootstrap samples
#' @param omega.freq.method A character string for the method of frequentist omega, either "pfa" or "cfa",
#' with "pfa" the interval is always bootstrapped
#' @param n.obs A number for the sample observations when a covariance matrix is supplied
#' and the factor model is calculated
#' @param alpha.int.analytic A logical for calculating the alpha confidence interval analytically
#' @param omega.int.analytic A logical for calculating the omega confidence interval analytically,
#' only works with cfa as the omega.freq.method
#' @param freq A logical for calculating the frequentist estimates
#' @param Bayes A logical for calculating the Bayesian estimates
#' @param para.boot A logical for calculating the parametric bootstrap, the default is the non-parametric
#' @param item.dropped A logical for calculating the if-item-dropped statistics
#' @param missing A string specifying the way to handle missing data, 'listwise' is self-explanatory,
#' 'pairwise' in the Bayesian paradigm means sampling the missing values as additional parameters
#' from the joint conditional distribution, in the frequentist paradigm this means using the 'pairwise' covariance
#' matrix and the full information ML method for omega
#' @param callback step count for external use
#'
#' @examples
#' summary(strel(asrm, estimates = "lambda2", n.chains = 2, n.iter = 200, n.boot = 200))
#' summary(strel(asrm, estimates = "lambda2", item.dropped = TRUE, n.chains = 2,
#' n.iter = 100, n.boot = 200))
#'
#'
#' @references{
#'   \insertRef{murphy2007}{Bayesrel}
#'   \insertRef{lee2007}{Bayesrel}
#' }
#' @importFrom grDevices adjustcolor recordPlot
#' @importFrom graphics arrows axis legend lines par plot text title
#' @importFrom methods is
#' @importFrom stats cov cov2cor density ecdf qnorm quantile rchisq rgamma rnorm runif sd var approxfun integrate
#' @importFrom Rdpack reprompt
#'
#' @useDynLib Bayesrel, .registration=TRUE
#' @importFrom Rcpp evalCpp
#'
#' @export
strel <- function(data = NULL, estimates = c("alpha", "lambda2", "glb", "omega"),
                  cov.mat = NULL,
                  interval = .95,
                  n.iter = 1000, n.burnin = 50, thin = 1, n.chains = 3,
                  n.boot = 1000,
                  omega.freq.method = "cfa",
                  n.obs = NULL,
                  alpha.int.analytic = TRUE,
                  omega.int.analytic = TRUE,
                  freq = TRUE, Bayes = TRUE,
                  para.boot = FALSE,
                  item.dropped = FALSE,
                  missing = "pairwise",
                  callback = function(){}) {

  default <- c("alpha", "lambda2", "lambda4", "lambda6", "glb", "omega")
  # estimates <- match.arg(arg = estimates, several.ok = T)
  mat <- match(default, estimates)
  estimates <- estimates[mat]
  estimates <- estimates[!is.na(estimates)]
  p <- NULL
  sum_res <- list()
  sum_res$call <- match.call()

  pairwise <- FALSE
  use.cases <- "everything"
  if (any(is.na(data))) {
    if (missing == "listwise") {
      pos <- which(is.na(data), arr.ind = T)[, 1]
      data <- data[-pos, ]
      ncomp <- nrow(data)
      sum_res$complete <- ncomp
      use.cases <- "complete.obs"
    } else if (missing == "pairwise") {
      pairwise <- T
      sum_res$miss_pairwise <- T
      use.cases <- "pairwise.complete.obs"
    } else return(warning("missing values in data detected, please remove and run again"))
  }


  if (!is.null(cov.mat)){
    if (is.null(n.obs))
      return(warning("number of observations (n.obs) needs to be specified when entering a covariance matrix"))
    if (sum(cov.mat[lower.tri(cov.mat)] != t(cov.mat)[lower.tri(cov.mat)]) > 0)
      return(warning("input matrix is not symmetric"))
    if (!("matrix" %in% class(try(solve(cov.mat),silent=TRUE))))
      return(warning("Data covariance matrix is not invertible"))
    data <- MASS::mvrnorm(n.obs, rep(0, ncol(cov.mat)), cov.mat, empirical = TRUE)
    colnames(data) <- colnames(cov.mat)
  }

  if (!("matrix" %in% class(try(solve(cov(data, use = use.cases)), silent=TRUE))))
      return(warning("Data covariance matrix is not invertible"))

  data <- scale(data, scale = F) # needed for Bayes omega

  if (Bayes) {
    sum_res$Bayes <- gibbsFun(data, estimates, n.iter, n.burnin, thin, n.chains, interval, item.dropped, pairwise,
                              callback)
    sum_res$n.iter <- n.iter
    sum_res$n.burnin <- n.burnin
    sum_res$thin <- thin
    sum_res$n.chains <- n.chains
  }


  if(freq){

    if (para.boot){
      sum_res$freq <- freqFun_para(data, n.boot, estimates, interval, omega.freq.method, item.dropped,
                                   alpha.int.analytic, omega.int.analytic, pairwise, callback)
    } else {
      sum_res$freq <- freqFun_nonpara(data, n.boot, estimates, interval, omega.freq.method, item.dropped,
                                    alpha.int.analytic, omega.int.analytic, pairwise, callback)
    }

    if(alpha.int.analytic) {sum_res$alpha.interval = "analytic"}
    if(omega.int.analytic) {sum_res$omega.interval = "analytic"}

    sum_res$n.boot <- n.boot
    sum_res$para.boot <- para.boot
  }


  sum_res$estimates <- estimates
  sum_res$interval <- interval
  sum_res$data <- data


  class(sum_res) = 'strel'
  return(sum_res)
}
