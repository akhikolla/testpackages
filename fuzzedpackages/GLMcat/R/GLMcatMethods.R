#' Summary of models
#' @description \code{summary} method for GLMcat objects.
#' @param object a GLMcat model
#' @param ... additional arguments affecting the summary produced.
#' @rdname summary
#' @export
summary.glmcat <- function(object, ...) {
  coef <- object$coefficients
  se <- object$stderr
  tval <- coef / se

  object$coefficients <- cbind(
    "Estimate" = coef,
    "Std. Error" = se,
    "z value" = tval,
    "Pr(>|z|)" = 2 * pnorm(-abs(tval))
  )
  colnames(object$coefficients) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  sum_ma <- object$coefficients
  printCoefmat(object$coefficients, P.values = TRUE, has.Pvalue = TRUE, ...)
  # cf src/stats/R/lm.R and case with no weights and an intercept
  # f <- object$fitted.values
  # r <- object$residuals
  # mss <- sum((f - mean(f))^2)
  # mss <- if (object$intercept) sum((f - mean(f))^2) else sum(f^2)
  # rss <- sum(r^2)
  #
  # object$r.squared <- mss/(mss + rss)
  # df.int <- if (object$intercept) 1L else 0L
  # n <- length(f)
  # rdf <- object$df
  # object$adj.r.squared <- 1 - (1 - object$r.squared) * ((n - df.int)/rdf)
  # class(object) <- "summary"
  # object
  # return(sum_ma)
}

#' Model coefficients
#' @description Extract model coefficients from a glmcat object.
#' @param object a GLMcat model.
#' @param na.rm TRUE for NA coefficients to be removed, default is FALSE.
#' @param ...	other arguments.
#' @rdname coef
#' @export
#' @examples
#' data(DisturbedDreams)
#' mod1 <- GLMcat(
#'   formula = Level ~ Age,
#'   ref_category = "Very.severe",
#'   data = DisturbedDreams, distribution = "logistic"
#' )
#' coef(mod1)
coef.glmcat <- function(object, na.rm = FALSE, ...) {
  if (na.rm) {
    coefs <- object$coefficients
    coefs[!is.na(coefs)]
  }
  else {
    object$coefficients
  }
}

#' Number of observations in a glmcat model
#' @description Extract the number of observations from a GLMcat model.
#' @param object a GLMcat model.
#' @param ...	other arguments.
#' @rdname nobs_glmcat
#' @export
#' @examples
#' data(DisturbedDreams)
#' mod1 <- GLMcat(
#'   formula = Level ~ Age,
#'   categories_order = c("Not.severe", "Severe.1", "Severe.2", "Very.severe"),
#'   data = DisturbedDreams, distribution = "logistic"
#' )
#' nobs_glmcat(mod1)
nobs_glmcat <- function(object, ...) {
  return(object$nobs_glmcat)
}

#' LogLikelihood glmcat models
#' @description Extract LogLikelihood for GLMcat models.
#' @rdname logLik
#' @param object a GLMcat model.
#' @param ...	other arguments.
#' @export
#' @examples
#' data(DisturbedDreams)
#' mod1 <- GLMcat(
#'   formula = Level ~ Age,
#'   categories_order = c("Not.severe", "Severe.1", "Severe.2", "Very.severe"),
#'   data = DisturbedDreams, distribution = "logistic"
#' )
#' logLik(mod1)
logLik.glmcat <- function(object, ...) {
  structure(object$LogLikelihood,
    df = object$df, nobs_glmcat = object$nobs_glmcat,
    class = "logLik"
  )
}
