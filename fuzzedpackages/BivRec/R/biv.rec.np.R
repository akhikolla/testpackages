#' Deprecated: Use bivrecNP
#'
#' @description
#' Deprecated function from the previous version. Use \verb{bivrecNP}.
#'
#' @importFrom stats as.formula
#' @importFrom stats na.omit
#' @importFrom stats quantile
#'
#' @param formula A formula with six variables indicating the bivariate alternating gap time
#'
#' response on the left of the ~ operator and the covariates on the right.
#'
#' The six variables on the left must have the same length and be given as
#'
#' \verb{id + episode + xij + yij + d1 + d2 ~ 1}, where:
#'
#' \itemize{
#'  \item \verb{id}: Vector of subject's unique identifier (i).
#'  \item \verb{episode}: Vector indicating the observation or episode (j) for a subject (i). This will determine order of events for each subject.
#'  \item \verb{xij}: Vector with the lengths of time spent in event of Type I for individual i in episode j.
#'  \item \verb{yij}: Vector with the lengths of time spent in event of Type II for individual i in episode j.
#'  \item \verb{d1}: Vector of censoring indicator corresponding to Type I gap times (xij): = 1 for uncensored, and = 0 for censored gap times.
#'  \item \verb{d2}: Vector of censoring indicator corresponding to Type II gap times (yij): = 1 for uncensored, and = 0 for censored gap times.
#' }
#'
#' @param data A data frame that includes all the vectors listed in the formula.
#' @param ai See details.
#' @param u1 A vector or single number to be used for the estimation of joint cdf P(Type I gap times \eqn{\le} u1, Type II gap times \eqn{\le} u2) in the nonparametric method.
#' @param u2 A vector or single number to be used for the estimation of joint cdf P(Type I gap times \eqn{\le} u1, Type II gap times \eqn{\le} u2) in the nonparametric method.
#' @param conditional A logical value. If TRUE, this function will calculate the conditional cdf for the Type II gap time given an interval of the Type I gap time and the bootstrap standard error and confidence interval at the specified confidence level. Default is FALSE.
#' @param given.interval A vector c(v1, v2) that must be specified if conditional = TRUE. The vector indicates an interval for the Type I gap time to use for the estimation of the cdf of the Type II gap time given this interval.
#' If given.interval = c(v1, v2), the function calculates P(Type II gap times \eqn{\le} y | v1 \eqn{\le} Type I gap times \eqn{\le} v2). The given values v1 and v2 must be in the range of gap times in the estimated marginal survival.
#' @param CI The level for confidence intervals the joint cdf, marginal survival and conditional cdf. Must be between 0.50 and 0.99. Default is 0.95.
#'
#' @details
#' \verb{ai} indicates a real non-negative function of censoring times to be used as weights in the nonparametric method. This variable can take on values of 1 or 2 which indicate:
#' \itemize{
#' \item \verb{ai=1}: the weights are simply 1 for all subjects, \eqn{a(Ci) = 1} (default).
#' \item \verb{ai=2}: the weight for each subject is the subject's censoring time, \eqn{a(Ci) = Ci}.
#' }
#'
#' @return See \verb{bivrecNP}.
#'
#' @references
#' Huang CY, Wang MC. (2005). Nonparametric estimation of the bivariate recurrence time distribution. Biometrics, 61: 392-402.
#' \url{https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1541-0420.2005.00328.x}
#'
#' @export

biv.rec.np <- function(formula, data, CI, ai, u1, u2, conditional, given.interval){

  .Deprecated("bivrecNP")

  if (missing(ai)) {ai<-1}
  if (missing(conditional)) {conditional <- FALSE}
  if (missing(CI)) {CI <- 0.95}
  if (missing(u1)) {u1 <- round(seq(quantile(xij, probs = 0.4), max(xij), length.out=5))}
  if (missing(u2)) {u2 <- round(seq(quantile(yij, probs = 0.4), max(yij), length.out=4))}

  ### PULL INFORMATION FROM PARAMETERS TO SEND TO REFORMAT
  variables <- all.vars(formula[[2]])

  ####Ensure unique identifiers are numeric
  iden <- eval(parse(text = paste("data$", variables[1], sep="")))
  iden.u <- unique(iden)
  new.id <- NULL
  if (class(iden)!="num") {
    if (class(iden)!="int") {
      for (i in 1:length(iden.u)){
        for (j in 1:length(iden)) {
          if (iden[j] == iden.u[i]){
            new.id=c(new.id,i)
          }
        }
      }
      data$new.id <- new.id
    }
  }
  data <- data[,-which(colnames(data)==variables[1])]
  colnames(data)[ncol(data)] = variables[1]

  ####extract vectors/data needed to send to biv.rec.reformat
  names <- paste("data$", variables, sep="")
  identifier <- eval(parse(text = names[1]))
  episode <- eval(parse(text = names[2]))
  xij <- eval(parse(text = names[3]))
  yij <- eval(parse(text = names[4]))
  if (length(names)==6) {
    c_indicatorX <- eval(parse(text = names[5]))
    c_indicatorY <- eval(parse(text = names[6]))
  } else {
    stop("Error: Must supply 6 arguments in left hand side of formula.")
  }

  bivrec_resp <- with(data, bivrecSurv(identifier, episode, xij,
                                       yij, c_indicatorX, c_indicatorY))

  ans <- bivrecNP(response = bivrec_resp, ai=ai, u1=u1, u2=u2,
                  conditional=conditional, given.interval=given.interval,
                  level=CI)

  print("See bivrecNP to view and plot results.")

  return(ans)


}

