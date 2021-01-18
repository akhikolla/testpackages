#' Factor a polynomial
#'
#' Factor a polynomial
#'
#' @param mpoly a character parseable by [mp()], an
#'   \code{mpoly} object, or a pointer to a polynomial in M2
#' @param code return only the M2 code? (default: \code{FALSE})
#' @param ... ...
#' @return a named list with elements \code{factor} (an
#'   \code{mpolyList} object) and \code{power}, an integer
#'   vector
#' @name factor_poly
#' @examples
#'
#' \dontrun{ requires Macaulay2 be installed and an interactive session
#'
#' ##### basic usage
#' ########################################
#'
#' ring("x", "y", coefring = "QQ")
#' factor_poly("x^4 - y^4")
#'
#' # reference function
#' factor_poly.("x^4 - y^4")
#'
#'
#' ##### different inputs
#' ########################################
#'
#' # factor_poly accepts mpoly objects:
#' # remember you must create the ring first!
#' (p <- mp("x^4 - y^4"))
#' factor_poly.(p)
#' factor_poly(p)
#' mp("(x-y) (x+y) (x^2+y^2)")
#'
#'
#'
#' ##### other examples
#' ########################################
#'
#' ring("x","y", "z", coefring = "QQ")
#' (p <- mp("(x^2 - y) (x^2 + y) (x + y)^2 (x - z)^2"))
#' factor_poly.(p)
#' factor_poly(p)
#'
#' (p <- mp("(x-1)^3 (y-1)^3"))
#' factor_poly.(p)
#' factor_poly(p)
#'
#' }



#' @rdname factor_poly
#' @export
factor_poly <- function (mpoly, code = FALSE) {

  # run m2
  args <- as.list(match.call())[-1]
  eargs <- lapply(args, eval, envir = parent.frame())
  pointer <- do.call(factor_poly., eargs)
  if(code) return(invisible(pointer))

  # parse output
  parsed_out <- m2_parse(pointer)

  # reformat and return
  list(
    factor = structure(
      lapply(parsed_out, function(.) mp(.[[1]])),
      class = "mpolyList"
    ),
    power = vapply(parsed_out, `[[`, integer(1), 2)
  )

}





#' @rdname factor_poly
#' @export
factor_poly. <- function (mpoly, code = FALSE, ...) {

  # prepare mpoly param
  if (is.m2_pointer(mpoly)) {
    mpoly_param <- m2_name(mpoly)
  } else if (is.mpoly(mpoly) || is.mpolyList(mpoly)) {
    mpoly_param <- mpolyList_to_m2_str(mpoly)
  } else {
    mpoly_param <- as.character(mpoly)
  }

  # create code
  m2_code <- sprintf("factor(%s)", mpoly_param)

  # message
  if (code) { message(m2_code); return(invisible(m2_code)) }

  # run m2 and return pointer
  m2.(m2_code)

}
