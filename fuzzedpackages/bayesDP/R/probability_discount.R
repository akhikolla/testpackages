
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @title Bayesian Discount Prior: Comparison Between Current and Historical Data
#' @description \code{probability_discount} can be used to estimate the posterior
#'   probability of the comparison between historical and current data in the
#'   context of a clinical trial with normal (mean) data.
#'   \code{probability_discount} is not used internally but is given for
#'   educational purposes.
#' @param mu scalar. Mean of the current data.
#' @param sigma scalar. Standard deviation of the current data.
#' @param N scalar. Number of observations of the current data.
#' @param mu0 scalar. Mean of the historical data.
#' @param sigma0 scalar. Standard deviation of the historical data.
#' @param N0 scalar. Number of observations of the historical data.
#' @param number_mcmc scalar. Number of Monte Carlo simulations. Default is 10000.
#' @param method character. Analysis method. Default value "\code{fixed}" estimates
#'   the posterior probability and holds it fixed. Alternative method "\code{mc}"
#'   estimates the posterior probability for each Monte Carlo iteration.
#'   See the the \code{bdpnormal} vignette \cr
#'   \code{vignette("bdpnormal-vignette", package="bayesDP")} for more details.
#' @details
#'   This function is not used internally but is given for educational purposes.
#'   Given the inputs,  the output is the posterior probability of the comparison
#'   between current and historical data in the context of a clinical
#'   trial with normal (mean) data.
#'
#' @return \code{probability_discount} returns an object of class "probability_discount".
#'
#' An object of class \code{probability_discount} contains the following:
#' \describe{
#'  \item{\code{p_hat}}{
#'    scalar. The posterior probability of the comparison historical data weight. If
#'    \code{method="mc"}, a vector of posterior probabilities of length
#'    \code{number_mcmc} is returned.}
#' }
#'
#' @references
#' Haddad, T., Himes, A., Thompson, L., Irony, T., Nair, R. MDIC Computer
#'   Modeling and Simulation working group.(2017) Incorporation of stochastic
#'   engineering models as prior information in Bayesian medical device trials.
#'   \emph{Journal of Biopharmaceutical Statistics}, 1-15.
#'
#' @examples
#' probability_discount(mu  = 0,   sigma = 1, N  = 100,
#'                      mu0 = 0.1, sigma0 = 1, N0 = 100)
#'
#' @rdname probability_discount
#' @import methods
#' @importFrom stats sd density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov pnorm
#' @aliases probability_discount,ANY-method
#' @export probability_discount
probability_discount <- setClass("probability_discount")

setGeneric("probability_discount",
           function(mu          = NULL,
                    sigma       = NULL,
                    N           = NULL,
                    mu0         = NULL,
                    sigma0      = NULL,
                    N0          = NULL,
                    number_mcmc = 10000,
                    method      = "fixed"){
             standardGeneric("probability_discount")
           })

setMethod("probability_discount",
          signature(),
          function(mu          = NULL,
                   sigma       = NULL,
                   N           = NULL,
                   mu0         = NULL,
                   sigma0      = NULL,
                   N0          = NULL,
                   number_mcmc = 10000,
                   method      = "fixed"){

  if(!(method %in% c("fixed", "mc")))
    stop("method entered incorrectly. Must be one of 'fixed' or 'mc'.")


  ### Preposterior of current mu using flat prior
  posterior_flat_sigma2 <- 1/rgamma(number_mcmc, (N - 1)/2, ((N - 1) * sigma^2)/2)
  s                     <- (posterior_flat_sigma2/((N-1)+1))^0.5
  posterior_flat_mu     <- rnorm(number_mcmc, mu, s)

  ### Posterior of historical data parameters using flat prior
  prior_sigma2 <- 1/rgamma(number_mcmc, (N0-1)/2, ((N0-1)*sigma0^2)/2)
  s0           <- (prior_sigma2/((N0-1)+1))^0.5
  prior_mu     <- rnorm(number_mcmc, mu0, s0)

  ### Test of model vs real
  p_hat <- mean(posterior_flat_mu < prior_mu)

  ### Transform probability to an intuitive value
  p_hat <- 2*ifelse(p_hat > 0.5, 1 - p_hat, p_hat)

  if(method == "mc"){
    Z     <- abs(posterior_flat_mu - prior_mu) / (s^2+s0^2)
    p_hat <- 2*(1-pnorm(Z))
  }

  return(p_hat)
})
