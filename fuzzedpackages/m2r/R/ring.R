#' Create a new ring in Macaulay2
#'
#' Create a new ring in Macaulay2
#'
#' @param vars vector of variable names
#' @param coefring coefficient ring (default: \code{"CC"})
#' @param order a term order (default: \code{"grevlex"})
#' @param code return only the M2 code? (default: \code{FALSE})
#' @param x formal argument for print method
#' @param ... ...
#' @return a reference to a Macaulay2 ring
#' @name ring
#' @examples
#'
#' \dontrun{ requires Macaulay2
#'
#' ##### basic usage
#' ########################################
#'
#' ring("x", "y")
#' ring("x", "y", coefring = "QQ")
#'
#'
#' ##### standard evaluation
#' ########################################
#'
#' ring_(c("x", "y"))
#' ring_(c("x", "y"), code = TRUE)
#'
#' (myring <- ring_(c("x1","x2","x3","y"), coefring = "QQ", order = "lex"))
#'
#' m2_name(myring)
#' m2_meta(myring, "vars")
#' m2_meta(myring, "coefring")
#' m2_meta(myring, "order")
#'
#' ##### other options
#' ########################################
#'
#' ring_.(c("x", "y"))
#' ring_.(c("x", "y"), code = TRUE)
#'
#' }









#' @rdname ring
#' @export
ring <- function(..., coefring = m2_coefrings(), order = m2_termorders(),
  code = FALSE
) {

  # grab args
  x <- list(vars = lapply(dots(...), eval, envir = parent.frame()))
  otherArgs <- as.list(match.call(expand.dots = FALSE))[-c(1:2)]

  # eval
  args <- lapply(c(x, otherArgs), eval)

  # run standard evaluation gb
  do.call("ring_", args)

}





#' @rdname ring
#' @export
ring. <- function(..., coefring = m2_coefrings(), order = m2_termorders(),
  code = FALSE
) {

  # grab args
  x <- list(vars = lapply(dots(...), eval, envir = parent.frame()))
  otherArgs <- as.list(match.call(expand.dots = FALSE))[-c(1:2)]

  # eval
  args <- lapply(c(x, otherArgs), eval)

  # run standard evaluation gb
  do.call("ring_.", args)

}








#' @rdname ring
#' @export
ring_ <- function(vars, coefring = m2_coefrings(), order = m2_termorders(),
  code = FALSE, ...
) {

  # arg checking
  coefring <- match.arg(coefring)
  order <- match.arg(order)

  # run ideal.
  args <- as.list(match.call())[-1]
  eargs <- lapply(args, eval, envir = parent.frame())
  pointer <- do.call(ring_., eargs)
  if(code) return(invisible(pointer))

  # run ring.
  # pointer <- ring_.(vars, coefring, order, code)
  # if(code) return(invisible(pointer))

  # construct R-side ring, class and return
  m2_structure(
    m2_name = m2_name(pointer),
    m2_class = "m2_polynomialring",
    m2_meta = list(
      vars = vars,
      coefring = coefring,
      order = order
    )
  )

}



#' @rdname ring
#' @export
ring_. <- function(vars, coefring = m2_coefrings(), order = m2_termorders(),
  code = FALSE, ...
) {

  # check args
  coefring <- match.arg(coefring)
  order <- match.arg(order)

  # caps order for M2
  order <- switch(order, lex = "Lex", glex = "GLex", grevlex = "GRevLex")

  # set glex order
  if (order == "GLex") {
    # {Weights => {1,1,1,1,1}, Lex => 4}
    m2order <- sprintf(
      "{Weights => {%s}, Lex => %d}",
      paste(rep("1",length(vars)), collapse = ","),
      length(vars)-1
    )
  } else {
    # {Lex => 5}
    m2order <- paste0("{", order, " => ", length(vars), "}")
  }

  # make ring name
  ringname <- name_and_increment("ring", "m2_ring_count")

  # sortedvars <- vars(reorder(mp(paste(vars, collapse = " ")), order = order))
  sortedvars <- vars

  # construct code and message
  line <- sprintf(
    "%s = %s[%s,MonomialOrder=>%s]",
    ringname, coefring, paste(sortedvars, collapse = ","), m2order
  )
  if(code) { message(line); return(invisible(line)) }

  # run m2
  ret <- m2.(line)

  m2_name(ret) <- ringname
  ret
}





# m2 coefficient rings currently supported
coefring_as_ring <- function(coefring) {

  m2_structure(
    m2_name = coefring,
    m2_class = "m2_polynomialring",
    m2_meta = list(
      vars = NULL,
      coefring = coefring,
      order = "grevlex"
    )
  )

}



m2_parse_object_as_function.m2_polynomialring <- function(x, params) {

  monoid <- params[[c(1)]]
  vars <- c(m2_meta(x, "vars"))
  order <- "grevlex"

  for (i in 1:length(monoid)) {
    if (!is.m2_option(monoid[[i]])) {
      vars <- c(vars, unlist(monoid[[i]]))
    } else if (monoid[[c(i,1)]] == "MonomialOrder") {
      for (j in 1:length(monoid[[c(i,2)]])) {
        if (
          is.m2_option(monoid[[c(i,2,j)]]) &&
          monoid[[c(i,2,j,1)]] %in% c("Lex", "Weights", "GRevLex")
        ) {

          order <- monoid[[c(i,2,j,1)]]

          # extra checking for glex
          if (
            order == "Weights" &&
            all(unlist(lapply(monoid[[c(i,2,j,2)]], function(x) x==1)))
          ) {
            order <- "glex"
            break()
          }

          order <- switch(order, Lex = "lex", GRevLex = "grevlex")

        }
      }
    }
  }

  m2_structure(
    m2_name = "",
    m2_class = "m2_polynomialring",
    m2_meta = list(
      vars = vars,
      coefring = m2_meta(x, "coefring"),
      order = order
    )
  )

}




# m2 coefficient rings currently supported
#' @rdname ring
#' @export
m2_coefrings <- function() c("CC", "RR", "QQ", "ZZ")




# m2 term orders currently supported
#' @rdname ring
#' @export
m2_termorders <- function() c("grevlex", "lex", "glex")




m2_ring_class_names <- function() {
  c(
    "Ring","PolynomialRing","QuotientRing",
    "InexactFieldFamily","InexactField"
  )
}


#' @rdname ring
#' @export
print.m2_polynomialring <- function(x, ...){

  s <- sprintf(
    "M2 Ring: %s[%s], %s order",
    m2_meta(x, "coefring"),
    paste(m2_meta(x, "vars"), collapse = ","),
    m2_meta(x, "order")
  )
  cat(s, "\n", sep = "")

  invisible(x)
}
