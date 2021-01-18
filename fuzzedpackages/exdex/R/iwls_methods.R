# =========================== coef.iwls() ================================== #

#' Extract Model Coefficients from an \code{"iwls"} object
#'
#' \code{coef} method for class \code{c("iwls", "exdex")}.
#'
#' @param object and object of class \code{c("kaps", "exdex")} returned from
#'   \code{\link{iwls}}.
#' @param ... Further arguments.  None are used.
#' @return A numeric scalar: the estimate of the extremal index \eqn{\theta}.
#' @export
coef.iwls <- function(object, ...) {
  if (!inherits(object, "exdex")) {
    stop("use only with \"exdex\" objects")
  }
  return(object$theta)
}

# =========================== nobs.iwls() ================================== #

#' Extract the Number of Observations from an \code{"iwls"} object
#'
#' \code{nobs} method for class \code{c("iwls", "exdex")}.
#'
#' @param object and object of class \code{c("iwls", "exdex")} returned from
#'   \code{\link{iwls}}.
#' @param ... Further arguments.  None are used.
#' @return A numeric scalar: the number of inter-exceedance times used in the
#'   fit.
#' @export
nobs.iwls <- function(object, ...) {
  if (!inherits(object, "exdex")) {
    stop("use only with \"exdex\" objects")
  }
  return(object$n_gaps)
}

# ============================ print.iwls() ================================== #

#' Print method for an \code{"iwls"} object
#'
#' \code{print} method for class \code{c("iwls", "exdex")}.
#'
#' @param x an object of class \code{c("iwls", "exdex")}, a result of
#'   a call to \code{\link{iwls}}.
#' @param digits The argument \code{digits} to \code{\link{print.default}}.
#' @param ... Additional arguments.  None are used in this function.
#' @details Prints the original call to \code{\link{iwls}}
#'   and the estimate of the extremal index \eqn{\theta}.
#' @return The argument \code{x}, invisibly, as for all
#'   \code{\link[base]{print}} methods.
#' @seealso \code{\link{iwls}} for maximum likelihood estimation of the
#'   extremal index \eqn{\theta} using the \eqn{K}-gaps model.
#' @export
print.iwls <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  if (!inherits(x, "exdex")) {
    stop("use only with \"exdex\" objects")
  }
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Convergence (0 means success):", x$conv, "\n\n")
  cat("Estimate of the extremal index theta:\n")
  print.default(format(coef(x), digits = digits), print.gap = 2L,
                quote = FALSE)
  return(invisible(x))
}
