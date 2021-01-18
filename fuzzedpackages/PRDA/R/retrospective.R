#############################
####    Retrospective    ####
#############################

#----    Retrospective    ----

#' Retrospective Design Analysis
#'
#' Given the hypothetical population effect size and the study sample size, the
#' function \code{retrospective()} performs a retrospective design analysis for
#' Pearson's correlation test between two variables or \emph{t}-test comparing
#' group means (Cohen's \emph{d}). According to the defined alternative
#' hypothesis and the significance level, inferential risks (i.e., Power level,
#' Type M error, and Type S error) are computed together with the critical
#' effect value (i.e., the minimum absolute effect size value that would result
#' significant).
#'
#'@param effect_size a numeric value or function (see Details) indicating the
#'  hypothetical population effect size.
#'@param sample_n1 a numeric value indicating the sample size of the first
#'  group.
#'@param sample_n2 a numeric value indicating the sample size of the second
#'  group. This argument is required when \code{test_method} is set to
#'  \code{"two_sample"} or \code{"welch"}. In the case of \code{test_method =
#'  "paired"}, set \code{sample_n2} equal to \code{sample_n1}. Whereas in the
#'  case of \code{test_method = "one_sample"}, set \code{sample_n2} to
#'  \code{NULL}. This argument is ignored for \code{test_method = "pearson"}.
#'  See Test methods section in Details.
#'@param test_method a character string specifying the test type, must be one of
#'  \code{"pearson"} (default, Pearson's correlation), \code{"two_sample"}
#'  (independent two-sample \emph{t}-test), \code{"welch"} (Welch's
#'  \emph{t}-test), \code{"paired"} (dependent \emph{t}-test for paired
#'  samples), or \code{"one_sample"} (one-sample \emph{t}-test). You can specify
#'  just the initial letters.
#'@param alternative a character string specifying the alternative hypothesis,
#'  must be one of \code{"two_sided"} (default), \code{"greater"} or
#'  \code{"less"}. You can specify just the initial letter.
#'@param sig_level a numeric value indicating the significance level on which
#'  the alternative hypothesis is evaluated.
#'@param ratio_sd a numeric value indicating the ratio between the standard
#'  deviation in the first group and in the second group. This argument is
#'  needed in the case of Welch's \emph{t}-test.
#'@param B a numeric  value indicating the number of iterations. Increase the
#'  number of iterations to obtain more stable results.
#'@param tl optional value indicating the lower truncation point if
#'  \code{effect_size} is defined as a function.
#'@param tu optional value indicating the upper truncation point if
#'  \code{effect_size} is defined as a function.
#'@param B_effect a numeric  value indicating the number of sampled effects
#'  if \code{effect_size} is defined as a function. Increase the number to
#'  obtain more stable results.
#'@param display_message a logical variable indicating whether to display or not
#'  the progress bar. Not that this applies only when \code{effect_size} is
#'  defined as a function.

#'
#'@return A list with class "design_analysis" containing the following
#'  components:
#'    \item{design_analysis}{a character string indicating the type of design
#'    analysis: "retrospective".}
#'    \item{call_arguments}{a list with all the arguments passed to the
#'    function and the raw function call.}
#'    \item{effect_info}{a list with all the information regarding the
#'    considered hypothetical population effect size. The list includes:
#'    \code{effect_type} indicating the type of effect; \code{effect_function}
#'    indicating the function from which effect are sampled or the string
#'    "single_value" if a single value was provided; \code{effect_summary}
#'    summary of the sampled effects; \code{effect_samples} vector with the
#'    sampled effects (or unique value in the case of a single value). if
#'    relevant \code{tl} and \code{tu} specifying the lower upper truncation
#'    point respectively.}
#'    \item{test_info}{a list with all the information regarding the test
#'    performed. The list includes: \code{test_method} character sting
#'    indicating the test method (i.e., "pearson", "one_sample", "paired",
#'    "two_sample", or "welch"); sample size (\code{sample_n1} and if relevant
#'    \code{sample_n2}), alternative hypothesis (\code{alternative}),
#'    significance level (\code{sig_level})  and  degrees of freedom (\code{df})
#'    of the statistical test; \code{critical_effect} the minimum absolute
#'    effect value that would result significant. Note that
#'    \code{critical_effect} in the case of \code{alternative = "two_sided"} is
#'    the absolute value and both positive and negative values should be
#'    considered.}
#'    \item{retrospective_res}{a data frame with the results of the design
#'    analysis. Columns names are \code{power}, \code{typeM}, and \code{typeS}.}
#'
#'@details Conduct a retrospective design analysis to evaluate inferential risks
#'  according to study design. A general overview is provided in the
#'  \code{vignette("retrospective")}.
#'
#'  \strong{Population effect size}
#'
#'  The hypothetical population effect size (\code{effect_size}) can be set to a
#'  single value or a function that allows sampling values from a given
#'  distribution. The function has to be defined as \code{function(n)
#'  my_function(n, ...)}, with only one single argument \code{n} representing
#'  the number of sampled values (e.g., \code{function(n) rnorm(n, mean = 0, sd
#'  = 1)}; \code{function(n) sample(c(.1,.3,.5), n, replace = TRUE)}). This
#'  allows users to define hypothetical effect size distribution according to
#'  their needs.
#'
#'  Argument \code{B_effect} allows defining the number of sampled effects.
#'  Users can access sampled effects in the \code{effect_info} list included in
#'  the output to evaluate if the sample is representative of their
#'  specification. Increase the number to obtain more accurate results but it
#'  will require more computational time (default is 1000). To avoid long
#'  computational times, we suggest adjusting \code{B} when using a function to
#'  define the hypothetical population effect size.
#'
#'  Optional arguments \code{tl} and \code{tu} allow truncating the sampling
#'  distribution specifying the lower truncation point and upper  truncation
#'  point respectively. Note that if \code{effect_type = "correlation"},
#'  distribution is automatically truncated between -1 and 1.
#'
#'  \strong{Test methods}
#'
#'  The function \code{retrospective()} performs a retrospective design analysis
#'  considering correlations between two variables or comparisons between group
#'  means.
#'
#'  In the case of a correlation, only Pearson's correlation between two
#'  variables is available, whereas Kendall's \emph{tau} and Spearman's
#'  \emph{rho} are not implemented. The \code{test_method} argument has to be
#'  set to \code{"pearson"} (default) and the \code{effect_size} argument is
#'  used to define the hypothetical population effect size in terms of Pearson's
#'  correlation coefficient (\eqn{\rho}). The \code{sample_n2} argument is
#'  ignored.
#'
#'  In the case of a comparison between group means, the \code{effect_size}
#'  argument is used to define the hypothetical population effect size in terms
#'  of Cohen's \emph{d} and the available \emph{t}-tests are selected specifying
#'  the argument \code{test_method}. For independent two-sample \emph{t}-test,
#'  use \code{"two_sample"} and indicate the sample size of the second group
#'  (\code{sample_n2}). For Welch's \emph{t}-test, use \code{"welch"} and
#'  indicate and indicate the sample size of the second group (\code{sample_n2})
#'  and the ratio between the standard deviation in the first group and in the
#'  second group (\code{ratio_sd}). For dependent \emph{t}-test for paired
#'  samples, use \code{"paired"} (\code{sample_n1} and \code{sample_n2} have to
#'  be equal). For one-sample \emph{t}-test, use \code{"one_sample"}
#'  (\code{sample_n2} has to be \code{NULL}).
#'
#'  \strong{Study design}
#'
#'  Study design can be further defined according to statistical test
#'  directionality and required \eqn{\alpha}-level using the arguments
#'  \code{alternative} and \code{sig_level} respectively.
#'
#'
#' @examples
#'
#' # Pearson's correlation
#' retrospective(effect_size = .3, sample_n1 = 25, test_method = "pearson")
#'
#' # Two-sample t-test
#' retrospective(effect_size = .3, sample_n1 = 25, sample_n2 = 35,
#'               test_method = "two_sample")
#' # Welch t-test
#' retrospective(effect_size = .3, sample_n1 = 25, sample_n2 = 35,
#'               test_method = "welch", ratio_sd = 1.5)
#' # Paired t-test
#' retrospective(effect_size = .3, sample_n1 = 25, sample_n2 = 25,
#'               test_method = "paired")
#' # One-sample t-test
#' retrospective(effect_size = .3, sample_n1 = 25, sample_n2 = NULL,
#'               test_method = "one_sample")
#'

#'
#'
#'
#' \donttest{
#' # Define effect_size using functions (long computational times)
#' # Remember to adjust B
#' retrospective(effect_size = function(n) rnorm(n, .3, .1), sample_n1 = 25,
#'               test_method = "pearson", tl = .15, B = 1e3)
#' retrospective(effect_size = function(n) rnorm(n, .3, .1), sample_n1 = 25,
#'               test_method = "one_sample", tl = .2, tu = .4, B = 1e3)
#' }
#'
#'@references Altoè, G., Bertoldo, G., Zandonella Callegher, C., Toffalini, E.,
#'  Calcagnì, A., Finos, L., & Pastore, M. (2020). Enhancing Statistical
#'  Inference in Psychological Research via Prospective and Retrospective Design
#'  Analysis. Frontiers in Psychology, 10.
#'  \url{https://doi.org/10.3389/fpsyg.2019.02893}
#'
#'  Bertoldo, G., Altoè, G., & Zandonella Callegher, C. (2020).
#'  Designing Studies and Evaluating Research Results: Type M and Type S Errors
#'  for Pearson Correlation Coefficient. Retrieved from
#'  \url{https://psyarxiv.com/q9f86/}
#'
#'  Gelman, A., & Carlin, J. (2014). Beyond Power Calculations: Assessing Type S
#'  (Sign) and Type M (Magnitude) Errors. Perspectives on Psychological Science,
#'  9(6), 641–651. \url{https://doi.org/10.1177/1745691614551642}
#'
#'
#'@export

retrospective <- function(effect_size,
                          sample_n1,
                          sample_n2 = NULL,
                          test_method = c("pearson", "two_sample", "welch",
                                          "paired", "one_sample"),
                          alternative = c("two_sided","less","greater"),
                          sig_level = .05,
                          ratio_sd = 1,
                          B = 1e4,
                          tl = -Inf,
                          tu = Inf,
                          B_effect = 1e3,
                          display_message = TRUE){



  #----    Save call    ----

  # Match arguments
  alternative <- match.arg(alternative)
  test_method <- match.arg(test_method)

  # eval effect_type
  effect_type <- eval_effect_type(test_method)

  # Save call
  design_analysis <- "retrospective"
  arguments_raw <- sys.call()
  call_arguments <- as.list(match_call(default = TRUE))[-1]

  # eval possible errors
  do.call(eval_arguments_retrospective,
          c(call_arguments,
            effect_type = effect_type))

  # Check sample_n2 for correlation
  if(effect_type == "correlation"){
    if(!is.null(sample_n2)){
      call_arguments["sample_n2"] <- list(NULL)
      warning("If 'test_method = pearson', argument 'ratio_n' is set to NULL")
    }
    sample_n2 <- NULL
  }

  #----    Evaluate effect size    ----

  effect_info <- eval_effect_size(effect_type = effect_type,
                                  effect_size = effect_size,
                                  tl = tl,
                                  tu = tu,
                                  B_effect = B_effect)

  effect_target <- effect_info$effect_summary[["Mean"]]

  #----    Get test info    ----

  # Evaluate test test_method
  # (use t.test or cor.tes() to evaluate possible errors)
  do.call(eval_test_method, c(call_arguments,
                              effect_type = effect_type,
                              effect_target = effect_target))

  # Compute df and critical value
  crit_values <- do.call(compute_critical_effect,
                         c(call_arguments,
                           effect_type = effect_type))

  test_info <- c(test_method = test_method,
                 sample_n1 = sample_n1,
                 sample_n2 = list(sample_n2), # list() used to deal with NULL
                 alternative = alternative,
                 sig_level = sig_level,
                 crit_values)

  #----    Retrospective analysis    ----

  retrospective_res <- do.call(simulate_analysis,
                               c(call_arguments,
                                 effect_type = effect_type,
                                 effect_info["effect_samples"]))

  #----    save results    ----
  design_fit <- structure(list(design_analysis = design_analysis,
                               call_arguments = c(call_arguments,
                                                  arguments_raw = arguments_raw),
                               effect_info = c(effect_type = effect_type,
                                               effect_info),
                               test_info = test_info,
                               retrospective_res = retrospective_res),
                          class = c("design_analysis","list"))

  return(design_fit)

}


#-----

