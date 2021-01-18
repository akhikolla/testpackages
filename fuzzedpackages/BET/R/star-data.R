#' Coordinates of Brightest Stars in the Night Sky
#'
#'This data set collects the galactic coordinates of the 256 brightest stars in the night sky (Perryman et al. 1997).
#'We consider the longitude (\code{x}) and sine latitude (\code{y}) here.
#'
#' @docType data
#'
#' @usage data(star)
#'
#' @keywords datasets
#'
#'
#' @examples
#' data(star)
#' BETs(cbind(star$x.raw, star$y.raw))
"star"
