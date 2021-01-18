#' Factor an integer into primes
#'
#' Factor an integer into primes
#'
#' @param n an integer or a polynomial
#' @param code return only the M2 code? (default: \code{FALSE})
#' @param ... ...
#' @return a data frame with integer columns \code{prime} and
#'   \code{power} or \code{m2_pointer} referencing the factorization
#'   in M2.
#' @name factor_n
#' @examples
#'
#' \dontrun{ requires Macaulay2
#'
#' ##### basic usage
#' ########################################
#'
#' 2^2 * 3^7 * 5^2 # = 218700
#' factor_n(218700)
#' factor_n.(218700)
#'
#' (df <- factor_n(218700))
#' df$prime
#' df$power
#' str(df)
#'
#'
#' factor_n(218700, code = TRUE)
#'
#'
#' ##### other options
#' ########################################
#'
#' (integer_pointer <- m2.("218700"))
#' m2_name(integer_pointer)
#' factor_n(integer_pointer, code = TRUE)
#' factor_n(integer_pointer)
#'
#'
#'
#' factor_n(3234432540)
#' factor_n(323443254223453)
#' factor_n(rpois(1, 1e4))
#'
#'
#' ##### known issues
#' ########################################
#'
#' # R doesn't handle big ints well. note in the following
#' # the m2 code number is different than the supplied number
#' factor_n(32344325422364353453, code = TRUE)
#'
#' # this can be circumvented by passing a string instead
#' factor_n("32344325422364353453", code = TRUE)
#'
#' # but if the factors are large, R can't handle the parsing well
#' factor_n("32344325422364353453")
#'
#' # here's a workaround:
#' factor_pointer <- factor_n.("32344325422364353453")
#' m2_meta(factor_pointer, "ext_str")
#' extract_factors <- function(pointer) {
#'   require(stringr)
#'   str <- m2_meta(pointer, "ext_str")
#'   str <- str_sub(str, 19, -2)
#'   str <- str_extract_all(str, "\\{[0-9]+,[0-9]+\\}")[[1]]
#'   str <- str_sub(str, 2, -2)
#'   str <- str_split(str, ",")
#'   df <- as.data.frame(t(simplify2array(str)))
#'   names(df) <- c("prime", "power")
#'   df
#' }
#' (df <- extract_factors(factor_pointer))
#'
#'
#' # using gmp (currently broken)
#' # factor_n("32344325422364353453", gmp = TRUE)
#' m2("11 * 479 * 6138607975396537")
#' 11 * 479 * 6138607975396537
#'
#' }
#'




#' @rdname factor_n
#' @export
factor_n <- function (n, code = FALSE, ...) {

  # run m2
  args <- as.list(match.call())[-1]
  eargs <- lapply(args, eval, envir = parent.frame())
  pointer <- do.call(factor_n., eargs)
  if (code) return(invisible(pointer))

  # parse output
  parsed_out <- m2_parse(pointer)

  # reformat and return
  list(
    prime = vapply(parsed_out, `[[`, integer(1), 1),
    power = vapply(parsed_out, `[[`, integer(1), 2)
  )

}


#' @rdname factor_n
#' @export
factor_n. <- function (n, code = FALSE, ...) {

  # arg checking
  # compute factor with m2, e.g. "2^2*67*97*9433"
  if (is.m2_pointer(n)) {
    param <- m2_name(n)
  } else if (is.numeric(n)) {
    param <- as.character(format(n, scientific = FALSE)) # for big Z's
  } else {
    param <- n
  }

  # create code and message
  m2_code <- sprintf("factor(%s)", param)
  if(code) { message(m2_code); return(invisible(m2_code)) }

  # run m2 and return pointer
  m2.(m2_code)

}



