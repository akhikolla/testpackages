

#' @title Checking of unit-norm data
#'
#' @description Utility for normalizing data without unit norms.
#'
#' @param x observations claimed to have unit norms. Either a matrix of size
#' \code{c(nx, p)} or a vector of length \code{p}.
#' @param warnings whether to show warnings if the normalization of
#' \code{x} happened.
#' @return A curated version of \code{x} with unit-norm observations and
#' possible zeros excluded.
#' @author Eduardo García-Portugués, Davy Paindaveine, and Thomas Verdebout.
#' @examples
#' check_unit_norm(c(sqrt(2), sqrt(2), 0) / 2)
#' check_unit_norm(1:3, warnings = FALSE)
#' check_unit_norm(rbind(c(0, 0, 0), c(0, 0, 1), 1:3, c(NA, 0, 1)),
#'                 warnings = FALSE)
#' @keywords internal
#' @export
check_unit_norm <- function(x, warnings = TRUE) {

  # Tolerance
  tol <- sqrt(.Machine$double.eps)

  # Name of x for informative warnings
  name_x <- deparse(substitute(x))

  # Check unit norms
  norms <- rowSums(rbind(x * x))
  if (any(abs(norms - 1) > tol, na.rm = TRUE)) {

    # Warning for normalization
    if (warnings) {

      warning(paste(name_x,
                  "does not have unit-norm rows, normalized internally."))

    }

    # Exclude zero observations
    zeros <- which(norms < tol)
    if (length(zeros) > 0) {

      # Warning for zero observations
      if (warnings) {

        warning(paste("Observation(s)", paste(zeros, collapse = ", "), "of",
                    name_x, "are zero, excluded internally."))

      }

      # Remove zeros depending on vector or matrix input of x
      if (length(norms) == 1) {

        stop(paste0("No remaining data in ", name_x,
                    " after excluding zeros."))

      } else {

        # Exclude zeros
        x <- x[-zeros, , drop = FALSE]
        norms <- norms[-zeros]
        if (nrow(x) == 0) {

          stop(paste0("No remaining data in ", name_x,
                      " after excluding zeros."))

        }

      }

    }

    # Normalize (propagates NAs)
    x <- x / sqrt(norms)

  }
  return(x)

}
