############################
####    PRDA version    ####
############################

#----    Package documentation    ----

#' PRDA: Prospective and Retrospective Design Analysis.
#'
#' Given an hypothetical value of effect size, {PRDA} performs a prospective
#' or retrospective design analysis to evaluate the inferential risks (i.e.,
#' power, Type M error, and Type S error) related to the study design. See
#' \code{vignette("PRDA")} for a brief introduction to \emph{Design
#' Analysis}.
#'
#' PRDA package can be used for Pearson's correlation between two variables
#' or mean comparisons (i.e., one-sample, paired, two-sample, and Welch's
#' t-test) considering an hypothetical value of \eqn{\rho} or Cohen's \emph{d}
#' respectively. See \code{vignette("retrospective")} for more details.
#'
#' @section Functions:
#' In {PRDA} there are two main functions:
#' \itemize{
#' \item{\strong{\code{retrospective()}}}. Given the hypothetical population
#' effect size and the study sample size, the function \code{retrospective()}
#' performs a retrospective design analysis. According to the defined
#' alternative hypothesis and the significance level, the inferential risks
#' (i.e., Power level, Type M error, and Type S error) are computed together
#' with the critical effect value (i.e., the minimum absolute effect size value
#' that would result significant). To know more about function arguments and
#' examples see the function documentation
#' \code{\link[PRDA:retrospective]{?retrospective}} and
#' \code{vignette("retrospective")}.
#'
#' \item{\strong{\code{prospective()}}}. Given the hypothetical population
#' effect size and the required power level, the function \code{prospective()}
#' performs a prospective design analysis. According to the defined alternative
#' hypothesis and the significance level, the required sample size is computed
#' together with the associated Type M error, Type S error, and the critical
#' effect value (i.e., the minimum absolute effect size value that would
#' result significant).  To know more about function arguments and examples see
#' the function documentation \code{\link[PRDA:prospective]{?prospective}}
#' and \code{vignette("prospective")}.
#' }
#'
#' @section Hypothetical Effect Size:
#' The hypothetical population effect size can be defined as a single value
#' according to previous results in the literature or experts indications.
#' Alternatively, {PRDA} allows users to specify a distribution of plausible
#' values to account for their uncertainty about the hypothetical population
#' effect size.  To know how to specify the hypothetical effect size according
#' to a distribution and an example of application see
#' \code{vignette("retrospective")}.
#'
#' @references Altoè, G., Bertoldo, G., Zandonella Callegher, C., Toffalini, E.,
#'  Calcagnì, A., Finos, L., & Pastore, M. (2020). Enhancing Statistical
#'  Inference in Psychological Research via Prospective and Retrospective Design
#'  Analysis. Frontiers in Psychology, 10.
#'  \url{https://doi.org/10.3389/fpsyg.2019.02893}
#'
#'  Bertoldo, G., Altoè, G., & Zandonella Callegher, C. (2020, June 15).
#'  Designing Studies and Evaluating Research Results: Type M and Type S Errors
#'  for Pearson Correlation Coefficient. Retrieved from
#'  \url{https://psyarxiv.com/q9f86/}
#'
#'  Gelman, A., & Carlin, J. (2014). Beyond Power Calculations: Assessing Type S
#'  (Sign) and Type M (Magnitude) Errors. Perspectives on Psychological Science,
#'  9(6), 641–651. \url{https://doi.org/10.1177/1745691614551642}
#'
#'
#' @importFrom stats rnorm t.test cor.test qt pt sd var cor median
#' @importFrom MASS mvrnorm
#' @importFrom utils capture.output
#' @importFrom pbapply pblapply
#'
#' @docType package
#' @name PRDA
NULL


#----    use_rcpp and rcppArmadillo    ----

## usethis namespace: start
#' @useDynLib PRDA, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
## usethis namespace: end
NULL


#----
