#' Set Macaulay2 ring
#'
#' \code{use_ring()} sets the default referent ring on the Macaulay2
#' side using the \code{use} function.
#'
#' @param ring a \code{m2_ring} (see [ring()]),
#'   \code{m2_ring_pointer} (see [ring.()]), or a character
#'   string containing the name of a ring in Macaulay2
#' @name use_ring
#' @examples
#'
#' \dontrun{ requires Macaulay2
#'
#'
#' ##### basic usage
#' ########################################
#'
#' ring("x", coefring = "QQ")
#' factor_poly("x^4 + 1")
#'
#' QQtxyz <- ring("t","x","y","z", coefring = "QQ")
#' gb("t^4 - x", "t^3 - y", "t^2 - z")
#'
#' ring("x", "y", "z", "t", coefring = "QQ")
#' gb("t^4 - x", "t^3 - y", "t^2 - z")
#'
#' use_ring(QQtxyz)
#' gb("t^4 - x", "t^3 - y", "t^2 - z")
#'
#' }




#' @export
#' @rdname use_ring
use_ring <- function(ring) {

  # arg check
  stopifnot(is.m2_polynomialring(ring) || is.m2_polynomialring_pointer(ring) || is.character(ring))

  # construct code
  if (is.m2(ring)) {
    m2_code <- sprintf("use %s", m2_name(ring))
  } else {
    m2_code <- sprintf("use %s", ring)
  }

  # run m2
  out <- m2.(m2_code)

  # return
  invisible(out)

}
