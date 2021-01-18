#' Default Language for \code{cnum}
#'
#' Function to check the default language for \code{cnum} functions.
#'
#' This package supports Traditional Chinese and Simplified Chinese. The
#' language can be specified with the \code{lang} parameter in every function,
#' with \code{"tc"} for Traditional Chinese and \code{"sc"} for Simplified
#' Chinese. The default is \code{"tc"}, but this can be changed by setting
#' \code{\link[base]{options}(cnum.lang = "sc")}.
#'
#' @return The default language for \code{cnum} functions.
#'
#' @seealso \itemize{ \item \link[=conversion]{Functions for conversion} \item
#'   \link[=tools]{Functions for detetction and extraction} }
#'
#' @examples
#' \donttest{
#' # Set the default language to Simplified Chinese
#' options(cnum.lang = "sc")
#' }
#' default_cnum_lang()
#'
#' @export
#'
default_cnum_lang <- function() {
  val <- getOption("cnum.lang")
  if (is.null(val))
    val <- "tc"
  if (!val %in% return_langs() || is.na(val) || length(val) != 1L)
    stop("unsupported language `", val, "` set by options(\"cnum.lang\").")
  val
}
