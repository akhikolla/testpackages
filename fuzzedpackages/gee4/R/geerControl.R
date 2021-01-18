namedList <- function(...) {
  L <- list(...)
  snm <- sapply(substitute(list(...)), deparse)[-1]
  if (is.null(nm <- names(L))) nm <- snm
  if (any(nonames <- nm == "")) nm[nonames] <- snm[nonames]
  setNames(L,nm)
}

#' @title Control of GEE-MCD/WGEE-MCD Model Fitting
#'
#' @description Construct control structures for GEE-MCD/WGEE-MCD model
#' fitting
#'
#' @param trace whether or not the value of the objective function and the
#'   parameters should be print on every trace'th iteration.
#' @param profile whether or not parameters should be estimated sequentially
#'   using the idea of profile likelihood.
#' @param ignore.const.term whether or not the constant term should be
#'   considered when calculating quasi-likelihood and QIC.
#' @param errorMsg whether or not the error message should be print.
#' @param use.weights.vec whether or not a user-specified weights.vec should be
#'   used
#'
#' @export geerControl
geerControl <- function(trace = FALSE, profile = TRUE, ignore.const.term = TRUE,
                        errorMsg = FALSE, use.weights.vec = FALSE)
{
  structure(namedList(trace, profile, ignore.const.term, errorMsg,
                      use.weights.vec),
            class = "geerControl")
}
