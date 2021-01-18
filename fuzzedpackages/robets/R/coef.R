#' Coef robets model
#' 
#' @param object An object of class \code{robets}.
#' @param ... Other undocumented arguments.
#'
#' @examples 
#' model <- robets(nottem)
#' coef(model)
#' @export
coef.robets <- function(object, ...)
{
  object$par
}