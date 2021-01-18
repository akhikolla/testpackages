HSC_PC <- function(...) stop("HSC_PC is an abstract class and cannot be instantiated!")

classify <- function(x, stream, ...) UseMethod("classify")
classify.default <- function(x, stream, ...) {
  stop(gettextf("classify mathod is not implemented for class '%s'.",paste(class(x), collapse=", ")))
}

classify.HSC_PC <- function(x, stream, ...) {
  stop(gettextf("classify mathod is not implemented for class '%s'.",paste(class(x), collapse=", ")))
}