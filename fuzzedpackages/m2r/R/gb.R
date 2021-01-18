#' Compute a Grobner basis with Macaulay2
#'
#' Compute a Grobner basis with Macaulay2
#'
#' \code{gb} uses nonstandard evaluation; \code{gb_} is the standard evaluation
#' equivalent.
#'
#' @param x a character vector of polynomials to be parsed by [mp()], a
#'   \code{mpolyList} object, an [ideal()] or pointer to an ideal
#' @param control a list of options, see examples
#' @param raw_chars if \code{TRUE}, the character vector will not be parsed by
#'   [mp()], saving time (default: \code{FALSE}). the down-side is that the
#'   strings must be formated for M2 use directly, as opposed to for [mp()].
#'   (e.g. \code{"x*y+3"} instead of \code{"x y + 3"})
#' @param code return only the M2 code? (default: \code{FALSE})
#' @param ... ...
#' @return an \code{mpolyList} object of class \code{m2_grobner_basis} or a
#'   \code{m2_grobner_basis_pointer} pointing to the same. See [mpolyList()].
#' @seealso [mp()], [use_ring()]
#' @name gb
#' @examples
#'
#' \dontrun{ requires Macaulay2
#'
#'
#' ##### basic usage
#' ########################################
#'
#' # the last ring evaluated is the one used in the computation
#' ring("t","x","y","z", coefring = "QQ")
#' gb("t^4 - x", "t^3 - y", "t^2 - z")
#'
#' # here's the code it's running in M2
#' gb("t^4 - x", "t^3 - y", "t^2 - z", code = TRUE)
#'
#'
#'
#' ##### different versions of gb
#' ########################################
#'
#' # standard evaluation version
#' poly_chars <- c("t^4 - x", "t^3 - y", "t^2 - z")
#' gb_(poly_chars)
#'
#' # reference nonstandard evaluation version
#' gb.("t^4 - x", "t^3 - y", "t^2 - z")
#'
#' # reference standard evaluation version
#' gb_.(poly_chars)
#'
#'
#'
#' ##### different inputs to gb
#' ########################################
#'
#' # ideals can be passed to gb
#' I <- ideal("t^4 - x", "t^3 - y", "t^2 - z")
#' gb_(I)
#'
#' # note that gb() works here, too, since there is only one input
#' gb(I)
#'
#' # ideal pointers can be passed to gb
#' I. <- ideal.("t^4 - x", "t^3 - y", "t^2 - z")
#' gb_(I.)
#'
#' # setting raw_chars is a bit faster, because it doesn't use ideal()
#' gb("t^4 - x", "t^3 - y", "t^2 - z", raw_chars = TRUE, code = TRUE)
#' gb("t^4 - x", "t^3 - y", "t^2 - z", raw_chars = TRUE)
#'
#'
#'
#' ##### more advanced usage
#' ########################################
#'
#' # the control argument accepts a named list with additional
#' # options
#' gb_(
#'   c("t^4 - x", "t^3 - y", "t^2 - z"),
#'   control = list(StopWithMinimalGenerators = TRUE),
#'   code = TRUE
#' )
#'
#' gb_(
#'   c("t^4 - x", "t^3 - y", "t^2 - z"),
#'   control = list(StopWithMinimalGenerators = TRUE)
#' )
#'
#'
#'
#' ##### potential issues
#' ########################################
#'
#' # when specifying raw_chars, be sure to add asterisks
#' # between variables to create monomials; that's the M2 way
#' ring("x", "y", "z", coefring = "QQ")
#' gb("x y", "x z", "x", raw_chars = TRUE, code = TRUE) # errors without code = TRUE
#' gb("x*y", "x*z", "x", raw_chars = TRUE, code = TRUE) # correct way
#' gb("x*y", "x*z", "x", raw_chars = TRUE)
#'
#'
#'
#'
#'
#'
#'
#'
#'
#' }




#' @export
#' @rdname gb
gb <- function(..., control = list(), raw_chars = FALSE, code = FALSE) {

  # grab args
  x <- list(x = lapply(dots(...), eval, envir = parent.frame()))
  if(is.list(x) && (length(x) == 1) && is.m2_ideal(x[[c(1,1)]])) x <- x[[1]]
  if(is.list(x) && (length(x) == 1) && is.m2_ideal_pointer(x[[c(1,1)]])) x <- x[[1]]
  otherArgs <- as.list(match.call(expand.dots = FALSE))[-c(1:2)]

  # eval
  args <- lapply(c(x, otherArgs), eval)

  # run standard evaluation gb
  do.call("gb_", args)

}






#' @export
#' @rdname gb
gb. <- function(..., control = list(), raw_chars = FALSE, code = FALSE) {

  # grab args
  x <- list(x = lapply(dots(...), eval, envir = parent.frame()))
  if(is.list(x) && (length(x) == 1) && is.m2_ideal(x[[c(1,1)]])) x <- x[[1]]
  if(is.list(x) && (length(x) == 1) && is.m2_ideal_pointer(x[[c(1,1)]])) x <- x[[1]]
  otherArgs <- as.list(match.call(expand.dots = FALSE))[-c(1:2)]

  # eval
  args <- lapply(c(x, otherArgs), eval)

  # run standard evaluation gb
  do.call("gb_.", args)

}







# value version of f (standard user version)
#' @rdname gb
#' @export
gb_ <- function(x, control = list(), raw_chars = FALSE, code = FALSE, ...) {

  # run m2
  args <- as.list(match.call())[-1]
  eargs <- lapply(args, eval, envir = parent.frame())
  pointer <- do.call(gb_., eargs)
  if(code) return(invisible(pointer))

  # parse output
  parsed_out <- m2_parse(pointer)

  # more parsing
  out <- m2_structure(
    parsed_out[1,],
    m2_name = m2_name(pointer),
    m2_class = "m2_grobner_basis",
    m2_meta = list(ideal = m2_meta(pointer, "ideal")),
    base_class = "mpolyList"
  )

  # return
  out

}








# reference version of f (returns pointer to m2 object)
#' @rdname gb
#' @export
gb_. <- function(x, control = list(), raw_chars = FALSE, code = FALSE, ...) {

  # make gb name
  gb_name <- name_and_increment("gb", "m2_gb_count")

  # create m2_code params
  if (raw_chars) {

    # x <- vapply(x, str_replace_all, character(1), "[ ]+", "*")
    ideal_name <- paste0("ideal(", paste0(x, collapse = ", "), ")")

  } else {

    if (is.m2_ideal(x)) {
      ideal_name <- m2_name(x)
    } else if (is.m2_ideal_pointer(x)) {
      ideal_name <- m2_name(x)
    } else {
      x <- do.call(ideal_., list(x = x))
      ideal_name <- m2_name(x)
    }

  }

  # prepare control string
  if (length(control) > 0) {
    control_string <- paste0(", ",
      paste(names(control), r_to_m2(unlist(control)), sep = " => ", collapse = ", ")
    )
  } else {
    control_string <- ""
  }

  # construct m2_code and message
  m2_code <- sprintf("%1$s = gb(%2$s%3$s); gens %1$s", gb_name, ideal_name, control_string)
  if(code) { message(m2_code); return(invisible(m2_code)) }

  # run m2
  out <- m2.(m2_code)

  # change name and return pointer
  m2_name(out) <- gb_name
  if(!raw_chars) m2_meta(out) <- c(m2_meta(out), list(ideal = ideal_name))
  out

}












r_to_m2 <- function(x) {
  if (length(x) > 1) return(vapply(x, r_to_m2, character(1)))
  if(is.logical(x)) return(tolower(as.character(x)))
  if(is.numeric(x)) return(x)
  if(is.character(x)) return(x)
  stop("unexpected input.")
}
