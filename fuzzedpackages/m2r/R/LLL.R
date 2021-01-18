#' LLL algorithm
#'
#' Macaulay2's implementation of the LLL algorithm. This implementation is still
#' under development and is currently untested.
#'
#' @param mat a matrix (integer entries)
#' @param control additional arguments to pass to LLL; see examples
#' @param code return only the M2 code? (default: \code{FALSE})
#' @seealso [m2_matrix()]
#' @name LLL
#'
#' @return an object of class \code{m2_matrix}
#' @examples
#'
#' \dontrun{ requires Macaulay2
#'
#' ##### basic usage
#' ########################################
#'
#' # example 1
#' M <- matrix(c(
#'   1, 1, 1, 1,
#'   2, 0, 3, 4,
#'   1, 0, 0, 0,
#'   0, 1, 0, 0,
#'   0, 0, 1, 0,
#'   0, 0, 0, 1
#' ), nrow = 6, byrow = TRUE)
#'
#' LLL(M)
#'
#'
#'
#'
#' # example 2 (wikipedia)
#' M <- matrix(c(
#'   1, -1, 3,
#'   1,  0, 5,
#'   1,  2, 6
#' ), nrow = 3, byrow = TRUE)
#'
#' LLL(M)
#'
#'
#' ##### control
#' ########################################
#'
#' M <- matrix(c(
#'   1, 1, 1, 1,
#'   2, 0, 3, 4,
#'   1, 0, 0, 0,
#'   0, 1, 0, 0,
#'   0, 0, 1, 0,
#'   0, 0, 0, 1
#' ), nrow = 6, byrow = TRUE)
#'
#' LLL(M, code = TRUE)
#' LLL(M, control = list(Strategy = "NTL"), code = TRUE)
#' LLL(M, control = list(Strategy = c("BKZ", "RealFP")), code = TRUE)
#'
#' LLL(M)
#' LLL(M, control = list(Strategy = "NTL"))
#' LLL(M, control = list(Strategy = c("BKZ", "RealFP")))
#' LLL(M, control = list(Strategy = c("BKZ", "RealQP")))
#'
#'
#'
#' # method timings with microbenchmark.  note they are roughly the same
#' # for this example matrix
#' microbenchmark::microbenchmark(
#'   "NTL" = LLL(M, control = list(Strategy = "NTL")),
#'   "BKZ_RealFP" = LLL(M, control = list(Strategy = c("BKZ", "RealFP"))),
#'   "BKZ_RealQP" = LLL(M, control = list(Strategy = c("BKZ", "RealQP"))),
#'   "BKZ_RealRR" = LLL(M, control = list(Strategy = c("BKZ", "RealRR")))
#' )
#'
#'
#'
#' ##### additional examples
#' ########################################
#'
#' LLL.(M)
#' LLL(M, code = TRUE)
#'
#'
#'
#' }
#'






#' @rdname LLL
#' @export
LLL <- function (mat, control = list(), code = FALSE) {

  # run m2
  args <- as.list(match.call())[-1]
  eargs <- lapply(args, eval, envir = parent.frame())
  pointer <- do.call(LLL., eargs)
  if(code) return(invisible(pointer))

  # parse output
  m2_parse(pointer)
}





#' @rdname LLL
#' @export
LLL. <- function (mat, control = list(), code = FALSE) {

  # arg checking
  # if (is.m2_matrix(mat)) mat <- mat$rmatrix
  if (is.m2_pointer(mat)) {
    param <- m2_name(mat)
  } else {
    if (!is.integer(mat)) stopifnot(all(mat == as.integer(mat)))
    param <- paste0("matrix", listify_mat(mat))
  }

  # prepare control string
  if (length(control) > 0) {
    if(length(names(control)) != 1 && names(control) != "Strategy") {
      stop("LLL only accepts control argument \"Strategy\".", call. = FALSE)
    }
    strategy <- listify(control$Strategy)
    if(length(control$Strategy) == 1L) strategy <- str_sub(strategy, 2, -2)
    control_string <- sprintf(", Strategy => %s", strategy)
  } else {
    control_string <- ""
  }

  # create code and message
  m2_code <- sprintf("LLL(%s%s)", param, control_string)
  if(code) { message(m2_code); return(invisible(m2_code)) }

  # run m2 and return pointer
  m2.(m2_code)
}




