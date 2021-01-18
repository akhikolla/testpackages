#' Construct a working test and implement an interim or the final analysis for a globally efficient adaptive design.
#'
#' This function is used as a routine by \code{\link{adaptive_analysis_norm_global}} and \code{\link{sample_size_norm_global}}.
#'
#' @param overall_sig_level Overall significance level in (0, 1). Default is 0.025.
#' @param work_beta Type II error probability of the working test. Default is 0.05.
#' @param cond_alpha Conditional Type I error probability in (0, 1). Default is 0.025.
#' @param cost_type_1_err The loss caused by erroneously rejecting the null hypothesis. If 0 is specified, the loss leading to the working test with its Type I error probability being \code{significance_level} will be calculated.  Default is 0.
#' @param cost_type_2_err The loss caused by erroneously accepting the null hypothesis. If 0 is specified, the loss will be set to the value of the maximum of the basic analysis schedule.  Default is 0.
#' @param prev_cost The value of \code{cost_type_1_err} of the working test in the analysis just before the current analysis.
#' @param min_effect_size The minimum effect size.  It should be positive. The working test will be constructed to have the power of \code{1 - work_beta} for this effect size.  Default is 1.
#' @param effect_size The effect size for which the probability of rejecting the null hypothesis will be calculated. If \code{cost_type_1_err = 0}, \code{effect_size} will be forced to be the null value, 0.  Default is 0.
#' @param basic_schedule_num The number of analysis of the working test. Default is 50.
#' @param basic_schedule_power Determine the intervals between analyses. Default is 2.
#' @param basic_schedule The basic analysis schedule arbitrarily specified by user.
#' @param prior_dist Prior distribution for effect sizes of \code{min_effect_size * 0:10 / 2}.
#' @param prev_time The time of the analysis just before the current analysis. Either \code{prev_time} or \code{next_time} should be 0. See the example below.
#' @param time The time of the current analysis.
#' @param next_time The time of the next analysis. Either \code{prev_time} or \code{next_time} should be 0. See the example below.
#' @param stat The value of the current test statistic. The value of \code{stat} should be 0 at \code{time = 0}.
#' @param input_check Indicate whether or not the arguments input by user contain invalid values.
#' @param out_process The values used in calculation will be output in addition to the main output. Default is FALSE.
#' @param simpson_div The precision determining the precision of numerical integration. The default value is 6.
#' @param tol_boundary The precision in calculation of the stopping boundary of the working test.
#' @param tol_cost The precision in calculation of the loss, \code{cost_type_1_error}.
#' @return List of values of the parameters specified, information of the working test, and the conditional probability of rejecting the null hypothesis.
#' @seealso
#' \code{\link{adaptive_analysis_norm_global}} for example of this function.
#' @importFrom stats dnorm pnorm qnorm
#' @export
work_test_norm_global <- function(
  overall_sig_level = 0.025,
  work_beta = 0.05,
  cond_alpha = overall_sig_level,
  cost_type_1_err = 0,
  cost_type_2_err = 0,
  prev_cost = 0,
  min_effect_size = 1,
  effect_size = 0,
  basic_schedule_num = 50,
  basic_schedule_power = 2,
  basic_schedule = 0,
  prior_dist = 0,
  prev_time = 0,
  time = 0,
  next_time = 0,
  stat = 0,
  input_check = TRUE,
  out_process = FALSE,
  simpson_div = 6,
  tol_boundary = 1e-8,
  tol_cost = 1e-8
  ) {

  if ( length(input_check) > 1 ) stop("'input_check' should be scalar.")
  if ( input_check ) {
    if ( length(overall_sig_level) != 1 ) stop("'overall_sig_level' should be scalar.")
    if ( length(work_beta) != 1 ) stop("'work_beta' should be scalar.")
    if ( length(cond_alpha) != 1 ) stop("'cond_alpha' should be scalar.")
    if ( length(cost_type_1_err) != 1 ) stop("'cost_type_1_err' should be scalar.")
    if ( length(cost_type_2_err) != 1 ) stop("'cost_type_2_err' should be scalar.")
    if ( length(prev_cost) != 1 ) stop("'prev_cost' should be scalar.")
    if ( length(min_effect_size) != 1 ) stop("'min_effect_size' should be scalar.")
    if ( length(effect_size) != 1 ) stop("'effect_size' should be scalar.")
    if ( length(basic_schedule_num) != 1 ) stop("'basic_schedule_num' should be scalar.")
    if ( length(basic_schedule_power) != 1 ) stop("'basic_schedule_power' should be scalar.")

    if ( length(basic_schedule) == 0 ) basic_schedule <- 0
    basic_schedule <- sort(unique(basic_schedule))
    if ( (length(basic_schedule) > 1) && (length(basic_schedule) < 2) ) stop("When 'basic_schedule' is specified, its length should be greater than 1.")

    if ( (!(length(prior_dist) == 1 && prior_dist[1] == 0)) && (length(prior_dist) != 11) ) stop("The length of 'prior_dist' should be 11. Otherwise, set c(0).")

    if ( length(prev_time) != 1 ) stop("'prev_time' should be scalar.")
    if ( length(time) != 1 ) stop("'time' should be scalar.")
    if ( length(next_time) != 1 ) stop("'next_time' should be scalar.")
    if ( length(stat) != 1 ) stop("'stat' should be scalar.")
    if ( length(out_process) != 1 ) stop("'out_process' should be scalar.")
    if ( length(simpson_div) != 1 ) stop("'simpson_div' should be scalar.")
    if ( length(tol_boundary) != 1 ) stop("'tol_boundary' should be scalar.")
    if ( length(tol_cost) != 1 ) stop("'tol_cost' should be scalar.")
  }

  return( work_test_norm_c(
    overall_sig_level,
    work_beta,
    cond_alpha,
    cost_type_1_err,
    cost_type_2_err,
    prev_cost,
    min_effect_size,
    effect_size,
    basic_schedule_num,
    basic_schedule_power,
    basic_schedule,
    prior_dist,
    prev_time,
    time,
    next_time,
    stat,
    input_check,
    out_process,
    simpson_div,
    tol_boundary,
    tol_cost
  ) )
}


#' Calculate sample size or power for a globally efficient adaptive design.
#'
#' \code{sample_size_norm_global} calculates the power if the time of the final
#' analysis is given and otherwise the sample size.
#' The computed power for \code{effect_size} is an approximate lower bound.
#' Sample size is also calculated on the basis of the bound.
#'
#' @param initial_test Designate the initial working test generated by \code{\link{work_test_norm_global}} function.
#' @param sample_size If \code{TRUE}, the function will return the sample size required by the globally efficient adaptive design to have the power of \code{target_power}. If \code{FALSE}, the function will return the power when the final interim analysis and the final analysis are conducted at \code{time} and \code{final_time}, respectively.
#' @param effect_size The effect size, on the basis of which the power or sample size calculation will be performed. In globally efficient designs, any real value is allowed.
#' @param time The time of the current analysis.
#' @param target_power The power, on the basis of which the sample size calculation will be performed.
#' @param final_time The time of the final analysis.
#' @param tol_sample_size The precision in calculation of the sample size.
#' @param input_check Indicate whether or not the arguments input by user contain invalid values.
#' @return It returns the sample size (when \code{sample_size = TRUE}) or the power (when \code{sample_size = FALSE}).
#' @seealso
#' \code{\link{adaptive_analysis_norm_global}} for example of this function.
#' @export
sample_size_norm_global <- function(
  initial_test = 0,
  sample_size = TRUE,
  effect_size = 0,
  time = 0,
  target_power = 0.8,
  final_time = 0,
  tol_sample_size = 1e-5,
  input_check = TRUE
  ) {

  if ( length(input_check) > 1 ) stop("'input_check' should be scalar.")
  if ( input_check ) {
    if ( length(sample_size) != 1 ) stop("'sample_size' should be scalar.")
    if ( length(effect_size) != 1 ) stop("'effect_size' should be scalar.")
    if ( length(time) != 1 ) stop("'time' should be scalar.")
    if ( length(target_power) != 1 ) stop("'target_power' should be scalar.")
    if ( length(final_time) != 1 ) stop("'final_time' should be scalar.")
    if ( length(tol_sample_size) != 1 ) stop("'tol_sample_size' should be scalar.")
  }

  return( sample_size_norm_c(
    initial_test,
    sample_size,
    effect_size,
    time,
    target_power,
    final_time,
    tol_sample_size,
    input_check
  ) )
}


#' Analyze data according to a globally efficient adaptive design.
#'
#' \code{adaptive_analysis_norm_global} performs an globally efficient adaptive test, 
#' a Frequentist adaptive test with the specified significance level
#' with full flexibility.
#' Normality with known variance is assumed for the test statistic
#' (more accurately, the test statistic is assumed to follow Brownian motion.)
#' Null hypothesis is fixed at 0 without loss of generality.
#' Exact p-value, median unbiased estimate and confidence limits proposed by Gao et al. (2013) can also be calculated.
#' For detailed illustration, see \code{vignette("adpss_ex")}.
#'
#' @param initial_test Designate the initial working test generated by \code{work_test_norm_global} function.
#' @param times The sequence of times (sample size or information level) at which analyses were conducted.
#' @param stats The sequence of test statistics.
#' @param costs The sequence of loss required to construct working tests. Specification is optional. Partial specification is allowed, in which non-specification may be represented by \code{0}.
#' @param final_analysis If \code{TRUE}, the result input will be regarded as complete (no more data will be obtained) and the significance level will be exhausted. If \code{FALSE}, the current analysis will be regarded as an interim analysis and the significance level will be preserved.
#' @param estimate If \code{TRUE}, p-value, median unbiased estimator and upper and lower confidence limits will be calculated.
#' @param ci_coef The confidence coefficient. Default is 0.95.
#' @param tol_est The precision of the calculated results.
#' @param input_check Indicate whether or not the arguments input by user contain invalid values.
#' @return It returns whether or not the result was statistically significant, a p-value and an exact confidence limits.
#' @references
#' Kashiwabara, K., Matsuyama, Y. An efficient adaptive design approximating fixed sample size designs. In preparation.
#' Gao, P., Liu, L., Mehta, C. (2013) Exact inference for adaptive group sequential designs. Stat Med 32: 3991-4005.
#' @examples
#' # Construct an initial working test
#' # Note: cost_type_1_err will be automatically calculated when 0 is specified.
#' init_work_test <- work_test_norm_global(min_effect_size = -log(0.65), cost_type_1_err=1683.458)
#'
#' # Sample size calculation
#' sample_size_norm_global(
#'   initial_test = init_work_test,
#'   effect_size = 11.11 / 20.02, # needs not be MLE
#'   time = 20.02,
#'   target_power = 0.75,
#'   sample_size = TRUE
#'   )
#' @seealso
#' \code{\link{work_test_norm_global}} and \code{\link{sample_size_norm_global}}.
#' @export
adaptive_analysis_norm_global <- function(
  initial_test = 0,
  times = 0,
  stats = 0,
  costs = 0,
  final_analysis = TRUE,
  estimate = TRUE,
  ci_coef = 0.95,
  tol_est = 1e-8,
  input_check = TRUE
  ) {

  if ( length(input_check) > 1 ) stop("'input_check' should be scalar.")
  if ( input_check ) {
    if ( length(final_analysis) != 1 ) stop("'final_analysis' should be scalar.")
    if ( length(estimate) != 1 ) stop("'estimate' should be scalar.")
    if ( length(ci_coef) != 1 ) stop("'ci_coef' should be scalar.")
    if ( length(tol_est) != 1 ) stop("'tol_est' should be scalar.")
  }

  est <- exact_est_norm_c(
    initial_test,
    times,
    stats,
    costs,
    final_analysis,
    estimate,
    ci_coef,
    tol_est,
    input_check
  )

  return( est )
}




