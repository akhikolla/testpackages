#' bnfit to cpts
#'
#' Convert a \code{bn.fit} object (a list of cpts from the bnlearn package)
#' into a list of ordinary array-like cpts
#'
#' @param x A \code{bn.fit} object
#' @export
bnfit_to_cpts <- function(x) {
  cpts <- lapply(x, function(e) {
    ep <- e$prob
    if (length(dim(ep)) == 1L) {
      en <- structure(list(dimnames(ep)[[1]]), names = e$node)
      dimnames(ep) <- en
    }
    ep
  })
}
