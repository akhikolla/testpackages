# =========================== coef.kgaps() ================================== #

#' Extract Model Coefficients from a \code{"kgaps"} object
#'
#' \code{coef} method for class \code{c("kgaps", "exdex")}.
#'
#' @param object and object of class \code{c("kaps", "exdex")} returned from
#'   \code{\link{kgaps}}.
#' @param ... Further arguments.  None are used.
#' @return A numeric scalar: the estimate of the extremal index \eqn{\theta}.
#' @export
coef.kgaps <- function(object, ...) {
  if (!inherits(object, "exdex")) {
    stop("use only with \"exdex\" objects")
  }
  return(object$theta)
}

# =========================== vcov.kgaps() ================================== #

#' Calculate Variance-Covariance Matrix for a \code{"kgaps"} object
#'
#' \code{vcov} method for class \code{c("kgaps", "exdex")}.
#'
#' @param object and object of class \code{c("kgaps", "exdex")} returned from
#'   \code{\link{kgaps}}.
#' @param ... Further arguments.  None are used.
#' @return A 1 by 1 numeric matrix containing the estimated variance of the
#'   estimator.
#' @export
vcov.kgaps <- function(object, ...) {
  if (!inherits(object, "exdex")) {
    stop("use only with \"exdex\" objects")
  }
  return(object$se ^ 2)
}

# =========================== nobs.kgaps() ================================== #

#' Extract the Number of Observations from a \code{"kgaps"} object
#'
#' \code{nobs} method for class \code{c("kgaps", "exdex")}.
#'
#' @param object and object of class \code{c("kgaps", "exdex")} returned from
#'   \code{\link{kgaps}}.
#' @param ... Further arguments.  None are used.
#' @return A numeric scalar: the number of inter-exceedance times used in the
#'   fit. If \code{x$inc_cens = TRUE} then this includes 2 censored
#'   observations.
#' @export
nobs.kgaps <- function(object, ...) {
  if (!inherits(object, "exdex")) {
    stop("use only with \"exdex\" objects")
  }
  return(object$ss$n_kgaps)
}

# ============================ print.kgaps() ================================== #

#' Print method for a \code{"kgaps"} object
#'
#' \code{print} method for class \code{c("kgaps", "exdex")}.
#'
#' @param x an object of class \code{c("kgaps", "exdex")}, a result of
#'   a call to \code{\link{kgaps}}.
#' @param digits The argument \code{digits} to \code{\link{print.default}}.
#' @param ... Additional arguments.  None are used in this function.
#' @details Prints the original call to \code{\link{kgaps}}
#'   and the estimate of the extremal index \eqn{\theta}.
#' @return The argument \code{x}, invisibly, as for all
#'   \code{\link[base]{print}} methods.
#' @seealso \code{\link{kgaps}} for maximum likelihood estimation of the
#'   extremal index \eqn{\theta} using the \eqn{K}-gaps model.
#' @seealso \code{\link{confint.kgaps}}: \code{confint} method for
#'   class \code{"kgaps"}.
#' @export
print.kgaps <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  if (!inherits(x, "exdex")) {
    stop("use only with \"exdex\" objects")
  }
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Estimate of the extremal index theta:\n")
  print.default(format(coef(x), digits = digits), print.gap = 2L,
                quote = FALSE)
  return(invisible(x))
}

# =============================== summary.kgaps =============================== #

#' Summary method for a \code{"kgaps"} object
#'
#' \code{summary} method for class \code{"kgaps"}
#'
#' @param object an object of class "kgaps", a result of a call to
#'   \code{\link{kgaps}}.
#' @param digits An integer. Used for number formatting with
#'   \code{\link[base:Round]{signif}}.
#' @param ... Additional arguments.  None are used in this function.
#' @return Returns a list containing the list element \code{object$call}
#'   and a numeric matrix \code{summary} giving the Estimate of the extremal
#'   index \eqn{\theta} and the estimated standard error (Std. Error).
#' @seealso \code{\link{kgaps}} for estimation of the extremal index
#'   \eqn{\theta} using a semiparametric maxima method.
#' @seealso \code{\link{confint.kgaps}} for estimation of confidence intervals
#'   for \eqn{\theta}.
#' @section Examples:
#' See the examples in \code{\link{kgaps}}.
#' @export
summary.kgaps <- function(object, digits = max(3, getOption("digits") - 3L),
                        ...) {
  if (!inherits(object, "exdex")) {
    stop("use only with \"exdex\" objects")
  }
  res <- object["call"]
  theta <- signif(object$theta, digits = digits)
  se <- signif(object$se, digits = digits)
  res$matrix <- cbind(`Estimate` = theta, `Std. Error` = se)
  rownames(res$matrix) <- c("theta")
  class(res) <- "summary.kgaps"
  return(res)
}

# ============================ print.summary.kgaps ============================ #

#' Print method for objects of class \code{"summary.kgaps"}
#'
#' \code{print} method for an object \code{x} of class \code{"summary.kgaps"}.
#'
#' @param x An object of class "summary.pm", a result of a call to
#'   \code{\link{summary.kgaps}}.
#' @param ... Additional arguments passed on to \code{\link{print.default}}.
#' @return Prints the numeric matrix \code{x$summary} returned from
#' \code{\link{summary.kgaps}}.
#' @seealso \code{\link{kgaps}} for estimation of the extremal index
#'   \eqn{\theta} using a semiparametric maxima method.
#' @seealso \code{\link{confint.kgaps}} for estimation of confidence intervals
#'   for \eqn{\theta}.
#' @section Examples:
#' See the examples in \code{\link{kgaps}}.
#' @export
print.summary.kgaps <- function(x, ...) {
  if (!inherits(x, "summary.kgaps")) {
    stop("use only with \"summary.kgaps\" objects")
  }
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  print(x$matrix, ...)
  invisible(x)
}
