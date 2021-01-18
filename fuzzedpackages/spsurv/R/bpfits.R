#' Bernstein Polynomial Based Proportional Hazards Model
#'
#' @export
#' @description Fits the BPPH model to time-to-event data.
#' @param formula a Surv object with time to event observations, right censoring status and explanatory terms.
#' @param degree Bernstein polynomial degree.
#' @param data a data.frame object.
#' @param approach Bayesian or maximum likelihood estimation methods, default is approach = "mle".
#' @param ... further arguments passed to or from other methods
#' @return An object of class `spbp`.
#' @seealso \code{\link[spsurv]{spbp}}, \code{\link[spsurv]{bppo}} and \code{\link[spsurv]{bpaft}} for other BP based models.
#' @examples
#'
#' library("spsurv")
#' data("veteran")
#'
#' fit <- bpph(Surv(time, status) ~ karno + factor(celltype),
#' data = veteran)
#'
#' summary(fit)
#'
#' @importFrom rstan stan sampling optimizing
#' @importFrom survival Surv frailty


bpph <- function(formula, degree, data,  approach = c("mle", "bayes"), ...){
  fit <- spbp.default(formula = formula,
               degree = degree,
               data = data,
               model = "ph",
               approach = match.arg(approach),
               ...)
  fit$call$formula <- match.call()$formula
  fit$call$data <- match.call()$data
  return(fit)
}

#' Bernstein Polynomial Based Proportional Odds Model
#'
#' @export
#' @description Fits the BPPO model to time-to-event data.
#' @param formula a Surv object with time-to-event observations, right censoring status and explanatory terms.
#' @param degree Bernstein polynomial degree.
#' @param data a data.frame object.
#' @param approach Bayesian or maximum likelihood estimation methods, default is approach = "mle".
#' @param ... further arguments passed to or from other methods
#' @return An object of class `spbp`.
#' @seealso \code{\link[spsurv]{spbp}}, \code{\link[spsurv]{bpph}} and \code{\link[spsurv]{bpaft}} for other BP based models.
#' @examples
#'
#'library("spsurv")
#' data("veteran")
#'
#' fit <- bppo(Surv(time, status) ~ karno + celltype,
#' data = veteran)
#'
#' summary(fit)
#'
#' @importFrom rstan stan sampling optimizing
#' @importFrom survival Surv frailty

bppo <- function(formula, degree, data,  approach = c("mle", "bayes"), ...){
  fit <- spbp.default(formula = formula,
               degree = degree,
               data = data,
               model = "po",
               approach = match.arg(approach),
               ...)
  fit$call$formula <- match.call()$formula
  fit$call$data <- match.call()$data
  return(fit)
}

#' Bernstein Polynomial Based Accelerated Failure Time Model
#'
#' @export
#' @description Fits the BPAFT model to time-to-event data.
#' @param formula a Surv object with time to event observations, right censoring status and explanatory terms.
#' @param degree Bernstein polynomial degree.
#' @param data a data.frame object.
#' @param approach Bayesian or maximum likelihood estimation methods, default is approach = "mle".
#' @param ... further arguments passed to or from other methods
#' @return An object of class `spbp`.
#' @seealso \code{\link[spsurv]{spbp}}, \code{\link[spsurv]{bpph}} and \code{\link[spsurv]{bppo}} for other BP based models.
#' @examples
#'
#' library("spsurv")
#' data("veteran")
#'
#' fit <- bpaft(Surv(time, status) ~ karno + celltype,
#' data = veteran)
#'
#' summary(fit)
#'
#' @importFrom rstan stan sampling optimizing
#' @importFrom survival Surv frailty

bpaft <- function(formula, degree, data, approach = c("mle", "bayes"), ...){
  fit <- spbp.default(formula =  formula,
               degree = degree,
               data = data,
               model = "aft",
               approach = match.arg(approach),
               ...)
  fit$call$formula <- match.call()$formula
  fit$call$data <- match.call()$data
  return(fit)
}
