#' Fit gastric emptying curves with Stan
#'
#' @param d A data frame with columns
#' \itemize{
#'   \item \code{rec} Record descriptor as grouping variable, e.g. patient ID
#'   \item \code{minute} Time after meal or start of recording.
#'   \item \code{vol} Volume of meal or stomach
#'  }
#' @param model_name Name of predefined model in
#' \code{gastempt/exec}. Use \code{stan_model_names()} to get a list
#' of available models.
#' @param lkj LKJ prior for kappa/tempt correlation, only required
#' for model linexp_gastro_2b. Values from 1.5 (strong correlation) to 50
#' (almost independent) are useful. See
#' \url{http://www.psychstatistics.com/2014/12/27/d-lkj-priors/} for examples.
#' @param student_df Student-t degrees of freedom for residual error;
#' default 5. Use 3 for strong outliers; values above 10 are close to gaussian
#' residual distribution.
#' @param init_r for stan, default = 0.2; Stan's own default is 2, which
#' often results in stuck chains.
#' @param chains for stan; default = 4. For debugging, use 1.
#' @param ... Additional parameter passed to \code{sampling}
#'
#' @return A list of class stan_gastempt with elements \code{coef, fit, plot}
#' \itemize{
#'   \item \code{coef} is a data frame with columns:
#'     \itemize{
#'       \item \code{rec} Record descriptor, e.g. patient ID
#'       \item \code{v0} Initial volume at t=0
#'       \item \code{tempt} Emptying time constant
#'       \item \code{kappa} Parameter \code{kappa} for
#'             \code{model = linexp}
#'       \item \code{beta} Parameter \code{beta} for \code{model = powexp}
#'       \item \code{t50} Half-time of emptying
#'       \item \code{slope_t50} Slope in t50; typically in units of ml/minute
#'  On error, coef is NULL
#'    }
#'   \item \code{fit} Result of class `stanfit`
#'   \item \code{plot} A ggplot graph of data and prediction. Plot of raw data is
#'      returned even when convergence was not achieved.
#'  }
#' @useDynLib gastempt, .registration = TRUE
#' @examples
#' \dontrun{
#'   dd = simulate_gastempt(n_records = 6, seed = 471)
#'   d = dd$data
#'   ret = stan_gastempt(d)
#'   print(ret$coef)
#' }
#' @import rstan
#' @importFrom utils capture.output
#' @export
stan_gastempt = function(d, model_name = "linexp_gastro_2b", lkj = 2,
                         student_df = 5L, init_r = 0.2, chains = 4,  ...){
  assert_that(all(c("record", "minute","vol") %in% names(d)))
  . = NULL # keep notes quiet
  is_linexp = grepl("linexp", model_name)
#  rstan_options(auto_write = TRUE)
#  options(mc.cores = parallel::detectCores())
  # Integer index of records
  d$record_i =  as.integer(as.factor(d$record))

  data = list(
    prior_v0 = median(d$vol[d$minute < 10]), # only required for _1x
    n = nrow(d),
    n_record = max(d$record_i),
    lkj = lkj,
    student_df = as.integer(student_df),
    record = d$record_i,
    minute = d$minute,
    volume = d$vol)
  mod = stanmodels[[model_name]]
  capture.output({
    #fit = suppressWarnings(sampling(mod, data = data))
    fit = suppressWarnings(sampling(mod, data = data, init_r = init_r,
                                    chains = chains))
  })
  cf = summary(fit)$summary[,1]
  if (is_linexp){
    coef = tibble(
      record = unique(d$record),
      v0 = cf[grep("v0\\[", names(cf))],
      kappa = cf[grep("^kappa", names(cf))],
      tempt = cf[grep("^tempt", names(cf))]
    ) %>% t50
    attr(coef, "mu_kappa") = cf["mu_kappa"]
    attr(coef, "sigma_kappa") = cf["sigma_kappa"]
  } else {
    # Assume powexp
    coef = tibble(
      record = unique(d$record),
      v0 = cf[grep("v0\\[", names(cf))],
      beta = cf[grep("^beta", names(cf))],
      tempt = cf[grep("^tempt", names(cf))]
    ) %>% t50
    attr(coef, "mu_beta") = cf["mu_beta"]
    attr(coef, "sigma_beta") = cf["sigma_beta"]
  }
  # Common attributes
  attr(coef, "sigma") = cf["sigma"]
  attr(coef, "lp") = cf["lp__"]
  #
  # Compute plot
  plot = ggplot(d, aes(x = minute, y = vol)) + geom_point() +
    facet_wrap(~ record) +
    expand_limits(x = 0, y = 0) # force zeroes to be visible
  minute = seq(min(d$minute), max(d$minute), length.out = 51)
  if (is_linexp){
    title = paste0("Stan fitted linexp function, model ", model_name)
    pred_func = linexp
  } else {
    ### TODO Not tested
    title = paste0("Fitted powexp function, model ", model_name)
    pred_func = powexp
  }
  newdata  = coef %>%
    rowwise() %>%
    do({
      vol = pred_func(minute, pars = . )
      tibble(record = .$record, minute = minute, vol = vol)
    })

  plot = plot + geom_line(data = newdata, col = "#006400")  +
    ggtitle(title, subtitle = comment(d))

  # Assemble return
  ret = list(coef = coef, fit = fit, plot = plot)
  class(ret) = "stan_gastempt"
  ret
}

#' Extract coefficients from stan_gastempt result
#'
#' @param object Result of a call to stan_gastempt
#' @param ... other arguments
#'
#' @return a data frame with coefficients. See \code{\link{nlme_gastempt}} for an example.
#' @export
coef.stan_gastempt = function(object, ...){
  Call = match.call(expand.dots = TRUE)
  sigdig = as.integer(Call[["signif"]])
  cf = object$coef
  if (!is.null(Call[["signif"]])) {
    cf[,-1] = lapply(cf[,-1], signif, sigdig)
  }
  cf
}

#' Plot data points and fit curve of an stan_gastempt fit
#'
#' @param x Result of a call to stan_gastempt
#' @param ... other arguments
#'
#' @return a ggplot object. Use \code{print()} if used
#' non-interactively to show the curve
#' @method plot stan_gastempt
#' @export
plot.stan_gastempt = function(x, ...){
  x$plot
}

if (FALSE) {
  library(gastempt)
  library(rstan)
  library(dplyr)
  library(assertthat)
  dd = simulate_gastempt(n_records = 6, seed = 471)
  d = dd$data
  ret = stan_gastempt(d)
  coef = ret$coef
  plot(ret$plot)
}
