#' An S4 class to represent a Hyper parameter.
#'
#' @import methods
#' @slot alpha1 A numeric value
#' @slot alpha2 A numeric value
#' @slot delta A numeric value
#' @slot ggamma A numeric value
#' @slot bbeta A numeric value
#'
#' @examples
#' new("Hparam", alpha1 = 1, alpha2 = 2, bbeta = 3, delta = 4, ggamma = 5)
#' @export
setClass(
  "Hparam",
  slots = c(
    alpha1 = "numeric",
    alpha2 = "numeric",
    delta  = "numeric",
    ggamma = "numeric",
    bbeta  = "numeric"
  ),
  prototype = list(
    alpha1 = numeric(),
    alpha2 = numeric(),
    delta  = numeric(),
    ggamma = numeric(),
    bbeta  = numeric()
  )
)


#' @export
setValidity("Hparam", function(object) {
  if (object@alpha1 < 0 |
    object@alpha2 < 0 |
    object@delta < 0 |
    object@ggamma < 0 |
    object@bbeta < 0) {
    "Hyperparameter should be non-negative!"
  }
})
