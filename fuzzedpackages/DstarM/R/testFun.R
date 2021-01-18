#' Test fun.density with lower and upper bounds
#'
#' @param fun.density A density function to be evaluated.
#' @param args Additional arguments for fun.density.
#' @param lower Lower bounds of the parameter space with which \code{fun.density} can be evaluated.
#' @param upper Upper bounds of the parameter space with which \code{fun.density} can be evaluated.
#'
#' @return Returns TRUE if no errors occurred, otherwise returns an error message
#' @details A function that is called whenever a nondefault density function is passed to \code{DstarM}. It does some rough error checking.
#'
#' @export
#'
#' @examples
#' lower = c(.5, -6, .1, 0, 0)
#' upper = c(2, 6, .99, .99, 10)
#' args = list(t = seq(0, 5, .01), pars = lower, boundary = 'lower',
#' DstarM = TRUE)
#' testFun(fun.density = Voss.density, lower = lower, upper = upper,
#' args = args)
#' # TRUE

testFun <- function(fun.density, lower, upper, args = list()) {
  # veranderen naar iets met getPdf?
  test <- vector("list", 4)
  args$pars <- lower
  args$boundary <- "lower"
  test[[1]] <- try(do.call(fun.density, args), silent = TRUE)
  args$boundary <- "upper"
  test[[2]] <- try(do.call(fun.density, args), silent = TRUE)
  args$pars <- upper
  test[[3]] <- try(do.call(fun.density, args), silent = TRUE)
  args$boundary <- "lower"
  test[[4]] <- try(do.call(fun.density, args), silent = TRUE)
  if (!all(unlist(lapply(test, is.numeric)))) {
    stop("testing fun.density with lower and upper bounds resulted in non-numeric output.",
      call. = FALSE)
  } else if (any(unlist(lapply(test, anyNA)))) {
    stop("testing fun.density with lower and upper bounds resulted in NaN output.",
      call. = FALSE)
  } else if (any(lengths(test) != length(args$t))) {
    stop("testing fun.density with lower and upper bounds resulted in output of wrong length.",
      call. = FALSE)
  } else if (any(unlist(lapply(test, function(x) any(x < 0))))) {
    stop("testing fun.density with lower and upper bounds resulted in negative values.",
      call. = FALSE)
  } else {
    return(TRUE)
  }
}

