

#' Generate an lvec with (random) values
#'
#' @param n number of elements in result vector
#' @param fun function that generates values. Should accept a number of 
#'   elements to generate as its first argument. 
#' @param ... additional arguments are passed on to \code{fun}.
#' @param chunk_size the size of the chunks of values with which to fill the
#'   resulting \code{\link{lvec}}. When not given it uses the value of the 
#'   option 'chunk_size' (see \code{\link{options}}) otherwise the default 
#'   value.
#'
#' @return
#' Returns an \code{\link{lvec}} with length \code{n}. The type is determined 
#' by the type of values returned by \code{fun}.
#'
#' @examples
#' # generate an lvec with random normally distributed values with sd of 10
#' x <- generate(2E6, rnorm, sd = 10)
#' # generate lvec with random letters; use sample; expects n as its second
#' # argument, but we work around that by explicitly naming first argument x
#' y <- generate(2E6, sample, x = letters, replace = TRUE)
#' 
#' @export
generate <- function(n, fun, ..., chunk_size = 5E6) {
  if (missing(chunk_size)) chunk_size <- getOption("chunk_size", chunk_size)
  stopifnot(chunk_size >= 1)
  nchunks <- ceiling(n / chunk_size)
  start <- 1
  result <- NULL
  for (i in seq_len(nchunks)) {
    stop <- min(n, start + chunk_size - 1)
    nchunk <- stop - start + 1
    sample <- fun(nchunk, ...)
    if (is.null(result)) {
      result <- as_lvec(sample)
      length(result) <- n
    } else {
      lset(result, range = c(start, stop), values = sample)
    }
    if ((i %% 100) == 0) gc()
    start <- stop + 1
  }
  result
}
