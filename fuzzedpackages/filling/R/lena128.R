#' lena image at size of \eqn{(128 \times 128)}
#'
#' \emph{Lena} is probably one of the most well-known example in image processing and computer vision.
#' Well, here is a brief introduction on \href{http://www.cs.cmu.edu/~chuck/lennapg/lenna.shtml}{the story of Lena}.
#'
#'
#' @docType data
#' @usage data(lena128)
#' @format matrix of size \eqn{(128\times 128)}
#' @keywords datasets
#' @references Gonzalez, Rafael C. and Woods, Richard E. (2017) \emph{Digital Image Processing} (4th ed.). ISBN 0133356728.
#'
#' @source \href{http://sipi.usc.edu/database/?volume=misc}{USC SIPI Image Database}
#'
#' @examples
#' data(lena128)
#' image(lena128, col=gray((0:100)/100), axes=FALSE, main="lena128")
"lena128"
