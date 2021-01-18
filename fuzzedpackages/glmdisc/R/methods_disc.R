methods::setMethod(
  "show",
  methods::signature(object = "glmdisc"),
  function(object) {
    methods::show(object@best.disc[[1]])
  }
)


print.glmdisc <- function(object) {
  print(object@best.disc[[1]])
}


summary.glmdisc <- function(object) {
  summary(object@best.disc[[1]])
}


simple_roc <- function(labels, scores) {
  labels <- labels[order(scores, decreasing = TRUE)]
  labels <- as.numeric(as.character(labels))
  data.frame(TPR = cumsum(labels) / sum(labels), FPR = cumsum(!labels) / sum(!labels), labels)
}



#' Wrapper for \code{\link[caret]{contr.ltfr}}
#'
#' Simply redirects the call to \code{\link[caret]{contr.ltfr}}, this is done
#' to avoid errors when the caret namespace is not attached.
#'
#' @param ... \code{\link[caret]{contr.ltfr}} parameters.
#' @export
#' @keywords internal
#' @examples
#' contr.ltfr(2)
contr.ltfr <- function(...) {
  caret::contr.ltfr(...)
}
