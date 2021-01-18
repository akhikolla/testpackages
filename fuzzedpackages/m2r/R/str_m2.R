#' Give the structure of a Macaulay2 ring
#'
#' Give the structure of a Macaulay2 ring
#'
#' @param object An \code{m2} object
#' @param ... ...
#' @return Invisible the object passed in.
#' @export
#' @examples
#'
#' \dontrun{ requires Macaulay2
#'
#' a <- m2("1")
#'
#' R <- ring(c("x1", "x2", "x3"))
#' str_m2(R)
#' str_m2.default(R)
#'
#' }
#'
str_m2 <- function (object, ...) UseMethod("str_m2")



#' @export
str_m2.default <- function(object, ...) {
  cat("M2 Object\n")
  cat(sprintf("    Type : %s\n", class(object)[1]))
  # cat(sprintf("  R Name : %s\n", deparse(substitute(object))))
  cat(sprintf(" M2 Name : %s\n", m2_name(object)))
}



#' @export
str_m2.m2_polynomialring <- function(object, ...){

  cat("M2 Object\n")
  cat(sprintf("     Type : %s\n", class(object)[1]))
  # cat(sprintf("  R Name : %s\n", deparse(substitute(object))))
  cat(sprintf("  M2 Name : %s\n", m2_name(object)))
  cat(sprintf(" Coefring : %s\n", m2_meta(object, "coefring")))
  cat(sprintf("     Vars : %s\n", paste(m2_meta(object, "vars"), collapse = ", ")))
  cat(sprintf("    Order : %s\n", m2_meta(object, "order")))

  invisible(object)
}
