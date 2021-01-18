#' Create a new ideal in Macaulay2
#'
#' Create a new ideal in Macaulay2
#'
#' @param x a listing of polynomials. several formats are accepted, see
#'   examples.
#' @param ring the referent ring in Macaulay2
#' @param raw_chars if \code{TRUE}, the character vector will not be parsed by
#'   [mp()], saving time (default: \code{FALSE}). the down-side is that the
#'   strings must be formated for M2 use directly, as opposed to for [mp()].
#'   (e.g. \code{"x*y+3"} instead of \code{"x y + 3"})
#' @param code return only the M2 code? (default: \code{FALSE})
#' @param ideal an ideal object of class \code{m2_ideal} or
#'   \code{m2_ideal_pointer}
#' @param e1,e2 ideals for arithmetic
#' @param I,J ideals or objects parsable into ideals
#' @param ... ...
#' @return a reference to a Macaulay2 ideal
#' @name ideal
#' @examples
#'
#' \dontrun{ requires Macaulay2
#'
#'
#' ##### basic usage
#' ########################################
#'
#' ring("x", "y", coefring = "QQ")
#' ideal("x + y", "x^2 + y^2")
#'
#'
#'
#' ##### different versions of gb
#' ########################################
#'
#' # standard evaluation version
#' poly_chars <- c("x + y", "x^2 + y^2")
#' ideal_(poly_chars)
#'
#' # reference nonstandard evaluation version
#' ideal.("x + y", "x^2 + y^2")
#'
#' # reference standard evaluation version
#' ideal_.(poly_chars)
#'
#'
#'
#' ##### different inputs to gb
#' ########################################
#'
#' ideal_(   c("x + y", "x^2 + y^2") )
#' ideal_(mp(c("x + y", "x^2 + y^2")))
#' ideal_(list("x + y", "x^2 + y^2") )
#'
#'
#'
#' ##### predicate functions
#' ########################################
#'
#' I  <- ideal ("x + y", "x^2 + y^2")
#' I. <- ideal.("x + y", "x^2 + y^2")
#' is.m2_ideal(I)
#' is.m2_ideal(I.)
#' is.m2_ideal_pointer(I)
#' is.m2_ideal_pointer(I.)
#'
#'
#'
#' ##### ideal radical
#' ########################################
#'
#' I <- ideal("(x^2 + 1)^2 y", "y + 1")
#' radical(I)
#' radical.(I)
#'
#'
#'
#' ##### ideal dimension
#' ########################################
#'
#' I <- ideal_(c("(x^2 + 1)^2 y", "y + 1"))
#' dimension(I)
#'
#' # dimension of a line
#' ring("x", "y", coefring = "QQ")
#' I <- ideal("y - (x+1)")
#' dimension(I)
#'
#' # dimension of a plane
#' ring("x", "y", "z", coefring = "QQ")
#' I <- ideal("z - (x+y+1)")
#' dimension(I)
#'
#'
#'
#' ##### ideal quotients and saturation
#' ########################################
#'
#' ring("x", "y", "z", coefring = "QQ")
#' (I <- ideal("x^2", "y^4", "z + 1"))
#' (J <- ideal("x^6"))
#'
#' quotient(I, J)
#' quotient.(I, J)
#'
#' saturate(I)
#' saturate.(I)
#' saturate(I, J)
#' saturate(I, mp("x"))
#' saturate(I, "x")
#'
#'
#' ring("x", "y", coefring = "QQ")
#' saturate(ideal("x y"), "x^2")
#'
#' # saturation removes parts of varieties
#' # solution over R is x = -1, 0, 1
#' ring("x", coefring = "QQ")
#' I <- ideal("(x-1) x (x+1)")
#' saturate(I, "x") # remove x = 0 from solution
#' ideal("(x-1) (x+1)")
#'
#'
#'
#' ##### primary decomposition
#' ########################################
#'
#' ring("x", "y", "z", coefring = "QQ")
#' I <- ideal("(x^2 + 1) (x^2 + 2)", "y + 1")
#' primary_decomposition(I)
#' primary_decomposition.(I)
#'
#' I <- ideal("x (x + 1)", "y")
#' primary_decomposition(I)
#'
#' # variety = z axis union x-y plane
#' (I <- ideal("x z", "y z"))
#' dimension(I) # =  max dimension of irreducible components
#' (Is <- primary_decomposition(I))
#' dimension(Is)
#'
#'
#'
#' ##### ideal arithmetic
#' ########################################
#'
#' ring("x", "y", "z", coefring = "RR")
#'
#' # sums (cox et al., 184)
#' (I <- ideal("x^2 + y"))
#' (J <- ideal("z"))
#' I + J
#'
#' # products (cox et al., 185)
#' (I <- ideal("x", "y"))
#' (J <- ideal("z"))
#' I * J
#'
#' # equality
#' (I <- ideal("x", "y"))
#' (J <- ideal("z"))
#' I == J
#' I == I
#'
#' # powers
#' (I <- ideal("x", "y"))
#' I^3
#'
#' }










#' @rdname ideal
#' @export
ideal <- function(..., raw_chars = FALSE, code = FALSE) {

  # grab args
  x <- list(x = lapply(dots(...), eval, envir = parent.frame()))
  otherArgs <- as.list(match.call(expand.dots = FALSE))[-c(1:2)]

  # eval
  args <- lapply(c(x, otherArgs), eval)

  # run standard evaluation gb
  do.call("ideal_", args)

}






#' @rdname ideal
#' @export
ideal. <- function(..., raw_chars = FALSE, code = FALSE) {

  # grab args
  x <- list(x = lapply(dots(...), eval, envir = parent.frame()))
  otherArgs <- as.list(match.call(expand.dots = FALSE))[-c(1:2)]

  # eval
  args <- lapply(c(x, otherArgs), eval)

  # run standard evaluation gb
  do.call("ideal_.", args)

}









#' @rdname ideal
#' @export
ideal_ <- function(x, raw_chars = FALSE, code = FALSE, ...) {

  # run ideal.
  args <- as.list(match.call())[-1]
  eargs <- lapply(args, eval, envir = parent.frame())
  pointer <- do.call(ideal_., eargs)
  if(code) return(invisible(pointer))

  # parse output
  parsed_out <- m2_parse(pointer)

  # construct R-side ideal, class and return
  m2_name(parsed_out) <- m2_name(pointer)
  parsed_out

}




#' @rdname ideal
#' @export
ideal_. <- function(x, raw_chars = FALSE, code = FALSE, ...) {

  # make ideal name
  ideal_name <- name_and_increment("ideal", "m2_ideal_count")

  # make ideal_param
  if (raw_chars) {
    ideal_param <- listify(x)
  } else {
    if (is.character(x)) {
      mpolys <- mpolyList_to_m2_str(mp(x))
      ideal_param <- listify(mpolys)
    } else if (is.list(x) && all(vapply(x, is.character, logical(1)))) {
      mpolys <- mpolyList_to_m2_str(mp(unlist(x)))
      ideal_param <- listify(mpolys)
    } else if (is.list(x) && all(vapply(x, is.mpoly, logical(1)))) {
      mpolys <- structure(x, class = "mpolyList")
      mpoly_strings_for_m2 <- mpolyList_to_m2_str(mpolys)
      ideal_param <- listify(mpoly_strings_for_m2)
    } else if (is.list(x) && all(vapply(x, is.numeric, logical(1)))) {
      # this is like c(mp("x y"), mp("x z"), mp("x"))
      stop(
        "you appear to have used c() on mpolys.\n",
        "  this input format is not accepted, use list() instead.",
        call. = FALSE
      )
    } else if (is.mpolyList(x)) {
      mpoly_strings_for_m2 <- mpolyList_to_m2_str(x)
      ideal_param <- listify(mpoly_strings_for_m2)
    } else {
      stop("unrecognized input x. see ?ideal", call. = FALSE)
    }
  }

  # construct m2_code and message
  m2_code <- sprintf("%s = ideal(%s)", ideal_name, ideal_param)
  if(code) { message(m2_code); return(invisible(m2_code)) }

  # run m2
  out <- m2.(m2_code)

  # change name and return
  m2_name(out) <- ideal_name
  out

}




m2_parse_function.m2_ideal <- function(x) {
  m2_structure(
    m2_name = "",
    m2_class = "m2_ideal",
    m2_meta = list(
      ring = m2_meta(x[[1]], "ring"),
      gens = structure(x[[1]][1,], class = "mpolyList")
    )
  )
}


m2_parse_function.m2_monomialideal <- function(x) {
  m2_structure(
    m2_name = "",
    m2_class = "m2_ideal",
    m2_meta = list(
      ring = m2_meta(x[[1]], "ring"),
      gens = structure(x[[1]][1,], class = "mpolyList")
    )
  )
}

# m2_parse_function.m2_ideal_list <- function(x) {
#   m2_structure(
#     m2_name = "",
#     m2_class = "m2_ideal_list",
#     m2_meta = list(
#       ring = m2_meta(x[[1]], "ring"),
#       gens = structure(x[[1]][1,], class = "mpolyList")
#     )
#   )
# }



#' @rdname ideal
#' @export
print.m2_ideal <- function(x, ...) {

  # from print.m2_polynomialring
  s <- sprintf(
    "ring %s[%s] (%s)",
    m2_meta(m2_meta(x, "ring"), "coefring"),
    paste(m2_meta(m2_meta(x, "ring"), "vars"), collapse = ","),
    m2_meta(m2_meta(x, "ring"), "order")
  )

  # ideal stuff
  if (length(m2_meta(x, "gens")) > 1) {
    with_gen <- "with generators :"
  } else {
    with_gen <- "with generator :"
  }
  cat("M2 Ideal of", s, with_gen, "\n")
  gens_strings <- print(m2_meta(x, "gens"), silent = TRUE)
  cat(paste("<", paste(gens_strings, collapse = ",  "), ">\n"))
  # cat(str_pad(gens_strings, nchar(gens_strings)+2, side = "left"), sep = "\n")
  invisible(x)

}






#' @rdname ideal
#' @export
print.m2_ideal_list <- function(x, ...) {

  # from print.m2_polynomialring
  s <- sprintf(
    "%s[%s] (%s)",
    m2_meta(m2_meta(x, "ring"), "coefring"),
    paste(m2_meta(m2_meta(x, "ring"), "vars"), collapse = ","),
    m2_meta(m2_meta(x, "ring"), "order")
  )

  cat("M2 List of ideals of", s, ":", "\n")
  lapply(x, function(ideal){
    gens_strings <- print(m2_meta(ideal, "gens"), silent = TRUE)
    cat(paste("<", paste(gens_strings, collapse = ",  "), ">\n"))
  })

  invisible(x)

}





#' @rdname ideal
#' @export
radical <- function(ideal, ring, code = FALSE, ...) {

  # run radical.
  args <- as.list(match.call())[-1]
  eargs <- lapply(args, eval, envir = parent.frame())
  pointer <- do.call(radical., eargs)
  if(code) return(invisible(pointer))

  # parse output
  parsed_out <- m2_parse(pointer)

  # construct R-side ideal, class and return
  m2_name(parsed_out) <- m2_name(pointer)
  m2_meta(parsed_out, "radical_of") <- ideal
  parsed_out

}



#' @rdname ideal
#' @export
radical. <- function(ideal, ring, code = FALSE, ...) {

  # arg check
  if (!is.m2_ideal(ideal) && !is.m2_ideal_pointer(ideal))
    stop("unrecognized input ideal. see ?ideal", call. = FALSE)

  # make radical ideal name
  radical_name <- name_and_increment("ideal", "m2_ideal_count")

  # check ring to be QQ or ZZ/p

  # construct code and message
  m2_code <- sprintf("%s = radical(%s)", radical_name, m2_name(ideal))
  if(!missing(ring)) m2_code <- paste(sprintf("use %s;", ring), m2_code)
  if(code) { message(m2_code); return(invisible(m2_code)) }

  # run m2
  out <- m2.(m2_code)

  # change name and return
  m2_name(out) <- radical_name
  out
}



#' @rdname ideal
#' @export
saturate <- function(I, J, code = FALSE, ...) {

  # run saturate.
  args <- as.list(match.call())[-1]
  eargs <- lapply(args, eval, envir = parent.frame())
  pointer <- do.call(saturate., eargs)
  if(code) return(invisible(pointer))

  # parse output
  parsed_out <- m2_parse(pointer)

  # construct R-side ideal, class and return
  m2_name(parsed_out) <- m2_name(pointer)
  m2_meta(parsed_out, "saturation_of") <- I
  parsed_out

}



#' @rdname ideal
#' @export
saturate. <- function(I, J, code = FALSE, ...) {

  # arg check
  if (!is.m2_ideal(I) && !is.m2_ideal_pointer(I))
    stop("unrecognized input ideal. see ?ideal", call. = FALSE)

  if (!missing(J)) {
    if (is.m2_ideal(J) || is.m2_ideal_pointer(J)) {
      second_param <- paste0(",", m2_name(J))
    } else if (is.mpoly(J)) {
      second_param <- paste0(",", mpolyList_to_m2_str(J))
    } else if (is.character(J)) {
      second_param <- paste0(",", mpolyList_to_m2_str(mp(J)))
    } else {
      stop("unrecognized input J. see ?ideal", call. = FALSE)
    }
  } else {
    second_param <- ""
  }

  # make saturation ideal name
  saturate_name <- name_and_increment("ideal", "m2_ideal_count")

  # construct code and message
  m2_code <- sprintf("%s = saturate(%s%s)", saturate_name, m2_name(I), second_param)
  if(code) { message(m2_code); return(invisible(m2_code)) }

  # run m2
  out <- m2.(m2_code)

  # change name and return
  m2_name(out) <- saturate_name
  out
}








#' @rdname ideal
#' @export
quotient <- function(I, J, code = FALSE, ...) {

  # run quotient.
  args <- as.list(match.call())[-1]
  eargs <- lapply(args, eval, envir = parent.frame())
  pointer <- do.call(quotient., eargs)
  if(code) return(invisible(pointer))

  # parse output
  parsed_out <- m2_parse(pointer)

  # construct R-side ideal, class and return
  m2_name(parsed_out) <- m2_name(pointer)
  m2_meta(parsed_out, "quotient_of") <- I
  parsed_out

}



#' @rdname ideal
#' @export
quotient. <- function(I, J, code = FALSE, ...) {

  # arg check
  if (!is.m2_ideal(I) && !is.m2_ideal_pointer(I))
    stop("unrecognized input ideal. see ?ideal", call. = FALSE)

  if (!missing(J)) {
    if (is.m2_ideal(J) || is.m2_ideal_pointer(J)) {
      second_param <- paste0(",", m2_name(J))
    } else if (is.mpoly(J)) {
      second_param <- paste0(",", mpolyList_to_m2_str(J))
    } else if (is.character(J)) {
      second_param <- paste0(",", mpolyList_to_m2_str(mp(J)))
    } else {
      stop("unrecognized input J. see ?ideal", call. = FALSE)
    }
  } else {
    second_param <- ""
  }

  # make saturation ideal name
  quotient_name <- name_and_increment("ideal", "m2_ideal_count")

  # construct code and message
  m2_code <- sprintf("%s = quotient(%s%s)", quotient_name, m2_name(I), second_param)
  if(code) { message(m2_code); return(invisible(m2_code)) }

  # run m2
  out <- m2.(m2_code)

  # change name and return
  m2_name(out) <- quotient_name
  out
}












#' @rdname ideal
#' @export
primary_decomposition <- function(ideal, code = FALSE, ...) {

  # run primary_decomposition.
  args <- as.list(match.call())[-1]
  eargs <- lapply(args, eval, envir = parent.frame())
  pointer <- do.call(primary_decomposition., eargs)
  if(code) return(invisible(pointer))

  # parse output
  parsed_out <- m2_parse(pointer)

  # construct m2 ideal list, class and return
  m2_name(parsed_out) <- m2_name(pointer)
  attr(parsed_out, "m2_meta") <- list(ring = m2_meta(ideal, "ring"))
  # m2_meta(parsed_out, "ring") <-  m2_meta(ideal, "ring")
  class(parsed_out) <- c("m2_ideal_list", class(parsed_out))
  parsed_out

}



#' @rdname ideal
#' @export
primary_decomposition. <- function(ideal, code = FALSE, ...) {

  # arg check
  if (!is.m2_ideal(ideal) && !is.m2_ideal_pointer(ideal))
    stop("unrecognized input ideal. see ?ideal", call. = FALSE)

  # make resulting structure name
  ideal_list_name <- name_and_increment("ideallist", "m2_ideal_list_count")

  # construct code and message
  m2_code <- sprintf("%s = primaryDecomposition(%s)", ideal_list_name, m2_name(ideal))
  if(code) { message(m2_code); return(invisible(m2_code)) }

  # run m2
  out <- m2.(m2_code)
  m2_name(out) <- ideal_list_name
  class(out) <- c("m2_ideal_list_pointer", class(out))

  # return
  out
}







#' @rdname ideal
#' @export
dimension <- function(ideal, code = FALSE, ...) {

  # arg check
  if (!(
    is.m2_ideal(ideal) || is.m2_ideal_pointer(ideal) ||
    is.m2_ideal_list(ideal) || is.m2_ideal_list_pointer(ideal)
  )) stop("unrecognized input ideal. see ?ideal", call. = FALSE)

  # construct code and message
  if(is.m2_ideal(ideal) || is.m2_ideal_pointer(ideal)) {
    m2_code <- sprintf("dim(%s)", m2_name(ideal))
  } else { # an ideal list or ideal list pointer
    m2_code <- sprintf("apply(%s, dim)", m2_name(ideal))
  }
  if(code) { message(m2_code); return(invisible(m2_code)) }

  # run m2
  out <- m2.(m2_code)

  # parse and return
  m2_parse(out)

}












#' @rdname ideal
#' @export
`+.m2_ideal` <- function(e1, e2) {

  # arg check
  if (!(
    (is.m2_ideal(e1) || is.m2_ideal_pointer(e1)) &&
    (is.m2_ideal(e2) || is.m2_ideal_pointer(e2))
  )) stop("unrecognized input ideal. see ?ideal", call. = FALSE)

  # construct code and message
  m2_code <- sprintf("%s + %s", m2_name(e1), m2_name(e2))

  # run m2
  out <- m2.(m2_code)

  # parse and return
  m2_parse(out)

}





#' @rdname ideal
#' @export
`*.m2_ideal` <- function(e1, e2) {

  # arg check
  if (!(
    (is.m2_ideal(e1) || is.m2_ideal_pointer(e1)) &&
    (is.m2_ideal(e2) || is.m2_ideal_pointer(e2))
  )) stop("unrecognized input ideal. see ?ideal", call. = FALSE)

  # construct code and message
  m2_code <- sprintf("%s * %s", m2_name(e1), m2_name(e2))

  # run m2
  out <- m2.(m2_code)

  # parse and return
  m2_parse(out)

}




#' @rdname ideal
#' @export
`==.m2_ideal` <- function(e1, e2) {

  # arg check
  if (!(
    (is.m2_ideal(e1) || is.m2_ideal_pointer(e1)) &&
    (is.m2_ideal(e2) || is.m2_ideal_pointer(e2))
  )) stop("unrecognized input ideal. see ?ideal", call. = FALSE)

  # construct code and message
  m2_code <- sprintf("%s == %s", m2_name(e1), m2_name(e2))

  # run m2
  out <- m2.(m2_code)

  # parse and return
  m2_parse(out)

}




#' @rdname ideal
#' @export
`^.m2_ideal` <- function(e1, e2) {

  # arg check
  if (!(
    (is.m2_ideal(e1) || is.m2_ideal_pointer(e1)) && is.numeric(e2)
  )) stop("unrecognized input ideal. see ?ideal", call. = FALSE)

  # construct code and message
  m2_code <- sprintf("%s^%s", m2_name(e1), e2)

  # run m2
  out <- m2.(m2_code)

  # parse and return
  m2_parse(out)

}
