#' Nonparametric Analysis of Bivariate Alternating Recurrent Event Gap Time Data
#'
#' @description
#' This function allows users to estimate the joint cumulative distribution function (cdf) for the two types of gap times (xij, yij), the marginal survival function for the Type I gap times (xij), and the conditional cdf for the Type II gap times (yij) given the Type I gap times (xij). See details for the estimation methods provided.
#'
#' @importFrom stats model.frame
#' @importFrom stats na.omit
#' @importFrom stats quantile
#' @importFrom stats model.matrix
#'
#' @param response A response object of class \verb{bivrecSurv}.
#' @param level The confidence level for confidence intervals for joint cdf, marginal survival probability and conditional cdf. Must be between 0.50 and 0.99. Default is 0.95.
#' @param ai See details.
#' @param u1 A vector or single number to be used for the estimation of joint cdf P(Type I gap times \eqn{\le} u1, Type II gap times \eqn{\le} u2) in the nonparametric method.
#' @param u2 A vector or single number to be used for the estimation of joint cdf P(Type I gap times \eqn{\le} u1, Type II gap times \eqn{\le} u2) in the nonparametric method.
#' @param conditional A logical value. If TRUE, this function will calculate the conditional cdf for the Type II gap time given an interval of the Type I gap time and the bootstrap standard error and confidence interval at the specified confidence level. Default is FALSE.
#' @param given.interval A vector c(v1, v2) that must be specified if conditional = TRUE. The vector indicates an interval for the Type I gap time to use for the estimation of the cdf of the Type II gap time given this interval.
#' If given.interval = c(v1, v2), the function calculates P(Type II gap times \eqn{\le} y | v1 \eqn{\le} Type I gap times \eqn{\le} v2). The given values v1 and v2 must be in the range of gap times in the estimated marginal survival.
#'
#' @details
#' \verb{ai} indicates a real non-negative function of censoring times to be used as weights in the nonparametric method. This variable can take on values of 1 or 2 which indicate:
#' \itemize{
#' \item \verb{ai=1} (default): the weights are simply 1 for all subjects, \eqn{a(Ci) = 1}.
#' \item \verb{ai=2}: the weight for each subject is the subject's censoring time, \eqn{a(Ci) = Ci}.
#' }
#'
#' Related methods: \verb{plot.bivrecNP}, \verb{head.bivrecNP}, \verb{print.bivrecNP}.
#'
#' @return A bivrecNP object that contains:
#' \itemize{
#'   \item \verb{joint_cdf}
#'   \item \verb{marginal_survival}
#'   \item \verb{conditional_cdf} (when conditional = TRUE)
#'   \item \verb{formula}
#'   \item \verb{ai}
#'   \item \verb{level}
#'   \item \verb{given.interval} (when conditional = TRUE)
#'   \item \verb{xij, yij}
#'   \item \verb{new_data}
#' }
#'
#' @references
#' Huang CY, Wang MC. (2005). Nonparametric estimation of the bivariate recurrence time distribution. Biometrics, 61: 392-402.
#' \url{https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1541-0420.2005.00328.x}
#'
#' @export
#' @examples
#'
#' \dontrun{
#'
#' library(BivRec)
#'
#' # Simulate bivariate alternating recurrent event data
#' set.seed(28)
#' sim_data <- simBivRec(nsize=100, beta1=c(0.5,0.5), beta2=c(0,-0.5),
#'             tau_c=63, set=1.1)
#' bivrecsurv_data <- with(sim_data, bivrecSurv(id, epi, xij, yij, d1, d2))
#' npresult <- bivrecNP(response = bivrecsurv_data, ai=1,
#'                      u1 = seq(2, 20, 1), u2 = seq(1, 15, 1), level=0.99)
#' head(npresult)
#' plot(npresult)
#'
#' #This is an example with longer runtime
#'
#'  npresult2 <- bivrecNP(response = bivrecsurv_data, ai=1,
#'                u1 = seq(2, 20, 1), u2 = seq(1, 15, 1), conditional = TRUE,
#'                given.interval = c(0, 10), level = 0.99)
#'  head(npresult2)
#'  plot(npresult2)
#' }
#'

bivrecNP <- function(response, ai, u1, u2, level, conditional, given.interval){

  x <- response

  if (!inherits(x, "bivrecSurv")) stop("Response must be a bivrecSurv object.")
  if (missing(ai)) {ai <- 1}
  if (missing(conditional)) {conditional <- FALSE}
  if (missing(level)) {CI <- 0.95} else {CI = level}

  if (CI > 0.99) {stop("Level is higher than 0.99")} else {
    if (CI<0.5) {stop("Level is less than 0.5")}
  }

  xij <- x$data4Creg$xij
  yij <- x$data4Creg$yij

  if (missing(u1)) {u1 <- round(seq(quantile(xij, probs = 0.4), quantile(xij, probs = 0.8), length.out=5))}
  if (missing(u2)) {u2 <- round(seq(quantile(yij, probs = 0.4), quantile(yij, probs = 0.8), length.out=4))}
  temp <- rep(u1, each = length(u2))
  temp2 <- rep(u2, length(u1))
  u <- cbind(u1=temp, u2=temp2)

  print("Estimating joint CDF and marginal survival")

  if (ai==1) {
    new_data = x$dat4np1
    forcdf <- new_data$forcdf
    formarg <- new_data$formarg
  } else {
    if (ai==2) {
      new_data = x$dat4np2
      forcdf <- new_data$forcdf
      formarg <- new_data$formarg
      } else {stop("ai must equal either 1 or 2.")}
    }

  cdf_res <- nonparam_cdf(forcdf, u, ai, CI)
  marg_res <- nonparam_marginal(formarg, CI)

  if (conditional==FALSE) {

    final_result <- list(joint_cdf = cdf_res, marginal_survival = marg_res, ai=ai,
                         xij=xij, yij=yij, new_data=new_data)
    final_result$level <- CI
    final_result$conditional <- conditional

    class(final_result) <- "bivrecNP"
    return(final_result)

  } else {

    if (missing(given.interval)) {
      print("Error: Missing given.interval while conditional=TRUE.")
      final_result <- list(joint_cdf = cdf_res, marginal_survival = marg_res, ai=ai)
      final_result$level <- CI
      final_result$conditional <- conditional

      class(final_result) <- "bivrecNP"
      return(final_result)

    } else {

      partial_result <- list(cdf = cdf_res, marginal_survival = marg_res,
                             ai=ai, new_data=new_data)

      ccdf_res <- nonparam_conditional(res=partial_result, given.interval, CI, yij)

      final_result <- list(joint_cdf = cdf_res, marginal_survival = marg_res,
                           conditional_cdf = ccdf_res, ai=ai, xij=xij, yij=yij, new_data=new_data)

      final_result$given.interval <- given.interval

      final_result$level <- CI
      final_result$conditional <- conditional

      class(final_result) <- "bivrecNP"
      return(final_result)
    }
  }
}
