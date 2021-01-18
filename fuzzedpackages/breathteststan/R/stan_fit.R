#' @title Bayesian Stan fit to 13C Breath Data
#' @description Fits exponential beta curves to 13C breath test series data using
#' Bayesian Stan methods. See
#' \url{https://menne-biomed.de/blog/breath-test-stan} for a comparision between
#' single curve, mixed-model population and Bayesian methods.
#'
#' @param data Data frame or tibble as created by \code{\link[breathtestcore]{cleanup_data}},
#' with mandatory columns \code{patient_id, group, minute} and \code{pdr}.
#' It is recommended to run all data through \code{\link[breathtestcore]{cleanup_data}} which
#' will insert dummy columns for \code{patient_id} and \code{minute} if the
#' data are distinct, and report an error if not. Since the Bayesian method
#' is stabilized by priors, it is possible to fit single curves.
#' @param dose Dose of acetate or octanoate. Currently, only one common dose
#' for all records is supported.
#' @param sample_minutes If mean sampling interval is < sampleMinutes, data are subsampled
#' using a spline algorithm
#' @param student_t_df When student_t_df < 10, the student distribution is used to
#' model the residuals. Recommended values to model typical outliers are from 3 to 6.
#' When student_t_df >= 10, the normal distribution is used.
#' @param chains Number of chains for Stan
#' @param iter Number of iterations for each Stan chain
#' @param model Name of model; use \code{names(stanmodels)} for other models.
#' @param seed Optional seed for rstan
#'
#'
#' @return A list of classes "breathteststanfit" and "breathtestfit" with elements
#' \itemize{
#'   \item {\code{coef} Estimated parameters as data frame in a key-value format with
#'    columns \code{patient_id, group, parameter, method} and \code{value}.
#'    Has an attribute AIC.}
#'    \item {\code{data}  The effectively analyzed data. If density of points
#'    is too high, e.g. with BreathId devices, data are subsampled before fitting.}
#'    \item {\code{stan_fit} The Stan fit for use with \code{shinystan::launch_shiny}
#'    or extraction of chains. }
#' }
#' @seealso Base methods \code{coef, plot, print}; methods from package
#'  \code{broom: tidy, augment}.
#' @examples
#' library(breathtestcore)
#' suppressPackageStartupMessages(library(dplyr))
#' d = breathtestcore::simulate_breathtest_data(n_records = 3) # default 3 records
#' data = breathtestcore::cleanup_data(d$data)
#' # Use more than 80 iterations and 4 chains for serious fits
#' fit = stan_fit(data, chains = 1, iter = 80)
#' plot(fit) # calls plot.breathtestfit
#' # Extract coefficients and compare these with those
#' # used to generate the data
#' options(digits = 2)
#' cf = coef(fit)
#' cf %>%
#'   filter(grepl("m|k|beta", parameter )) %>%
#'   select(-method, -group) %>%
#'   tidyr::spread(parameter, value) %>%
#'   inner_join(d$record, by = "patient_id") %>%
#'   select(patient_id, m_in = m.y, m_out = m.x,
#'          beta_in = beta.y, beta_out = beta.x,
#'          k_in = k.y, k_out = k.x)
#' # For a detailed analysis of the fit, use the shinystan library
#' \donttest{
#' library(shinystan)
#' # launch_shinystan(fit$stan_fit)
#' }
#' # The following plots are somewhat degenerate because
#' # of the few iterations in stan_fit
#' suppressPackageStartupMessages(library(rstan))
#' stan_plot(fit$stan_fit, pars = c("beta[1]","beta[2]","beta[3]"))
#' stan_plot(fit$stan_fit, pars = c("k[1]","k[2]","k[3]"))
#' stan_plot(fit$stan_fit, pars = c("m[1]","m[2]","m[3]"))
#'
#' @import rstan
#' @import rstantools
#' @import Rcpp
#' @import dplyr
#' @importFrom methods new
#' @useDynLib breathteststan, .registration = TRUE
#' @importFrom stats rnorm rlnorm
#' @importFrom utils capture.output
#' @importFrom stringr str_extract str_match
#' @importFrom tibble as_tibble
#' @importFrom purrr map_df
#' @importFrom stats na.omit quantile
#'
#' @export
#'
stan_fit = function(data, dose = 100, sample_minutes = 15, student_t_df = 10,
                    chains = 2, iter = 1000, model = "breath_test_1", seed = 4711) {

  # Avoid notes on CRAN
  value = pat_group = pat_group_i = NULL
  stat = estimate = . = k = key =  m = q_975 = NULL
  cm = comment(data)
  data = breathtestcore::subsample_data(data, sample_minutes)
  # Integer index of records
  data$pat_group_i =  as.integer(as.factor(data$pat_group))
  n_record = max(data$pat_group_i)

  data_list = list(
    n = nrow(data),
    n_record = n_record,
    dose = 100,
    student_t_df = student_t_df,
    pat_group_i = data$pat_group_i,
    minute = data$minute,
    pdr = data$pdr)

  # Note: as.array is required to handle the case of n_record = 1
  init = rep(list(list(
    m_raw = as.array(rnorm(n_record,0,.1)),
    mu_m = rnorm(1,40,2),
    sigma_m = abs(rnorm(1,2,.1)),

    k_raw = as.array(rnorm(n_record, 0,.1)),
    mu_k = rlnorm(1, -6,.1),
    sigma_k = abs(rnorm(1,0,.001)),

    beta_raw = as.array(rnorm(n_record, 0, .1)),
    mu_beta = rnorm(1, 2, 0.1),
    sigma_beta = abs(rnorm(1,.1,.1)),
    sigma = abs(rnorm(1,1,.1))
  )),chains)

  if (!exists("stanmodels"))
    stop("stanmodels not found")
  mod = stanmodels[[model]]
  if (is.null(mod))
    stop("Stan model", model,  "not found")
  options(mc.cores = min(chains, max(parallel::detectCores()/2, 1)))
  capture.output({fit = suppressWarnings(
    rstan::sampling(mod, data = data_list, init = init,
                    control = list(adapt_delta = 0.9),
                    seed = seed,
                    iter =  iter, chains = chains)
  )})

  # Extract required parameters
  cf = data.frame(pat_group_i = rep(1:n_record, each = chains*iter/2),
        m = as.vector(rstan::extract(fit, permuted = TRUE, pars = "m")$m),
        beta = as.vector(rstan::extract(fit, permuted = TRUE, pars = "beta")$beta),
        k = as.vector(rstan::extract(fit, permuted = TRUE, pars = "k")$k))
  # Compute derived quantities
  coef_chain = cf %>%
    mutate(
      t50_maes_ghoos = breathtestcore::t50_maes_ghoos(.),
      t50_maes_ghoos = breathtestcore::t50_maes_ghoos(.),
      tlag_maes_ghoos = breathtestcore::tlag_maes_ghoos(.),
      t50_maes_ghoos_scintigraphy = breathtestcore::t50_maes_ghoos_scintigraphy(.),
      t50_bluck_coward = breathtestcore::t50_bluck_coward(.),
      tlag_bluck_coward = breathtestcore::tlag_bluck_coward(.)
    ) %>%
    rename(m_exp_beta = m, k_exp_beta = k, beta_exp_beta = beta) %>%
    tidyr::gather(key, value, -pat_group_i) %>%
    na.omit() %>%
    ungroup()
  cf = coef_chain %>%
    group_by(pat_group_i, key) %>%
    summarize(
      estimate = mean(value),
      q_0275 = quantile(value, 0.0275),
      q_25 = quantile(value, 0.25),
      q_75 = quantile(value, 0.75),
      q_975 = quantile(value, 0.975)
    ) %>%
    ungroup() %>%
    left_join(unique(data[,c("pat_group_i", "pat_group", "patient_id", "group")]),
                by = "pat_group_i") %>%
    mutate(
      parameter = str_match(key, "k|m|beta|t50|tlag")[,1],
      method = str_match(key, "maes_ghoos_scintigraphy|maes_ghoos|bluck_coward|exp_beta")[,1]
    ) %>%
   select(-pat_group_i, -pat_group, -key)
  # Warning:
  # attributes are not identical across measure variables; they will be dropped
  cf = suppressWarnings(cf %>% tidyr::gather(stat, value, estimate:q_975))

  data = data %>% select(-pat_group, -pat_group_i) # only used locally
  ret = list(coef = cf, data = data, stan_fit = fit, coef_chain = coef_chain)
  comment(ret) = cm
  class(ret) = c("breathteststanfit", "breathtestfit")
  ret
}

