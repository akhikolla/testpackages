#' lena image at size of \eqn{(256 \times 256)}
#'
#' \emph{Lena} is probably one of the most well-known example in image processing and computer vision.
#' Well, here is a brief introduction on \href{http://www.cs.cmu.edu/~chuck/lennapg/lenna.shtml}{the story of Lena}.
#'
#'
#' @docType data
#' @usage data(lena256)
#' @format matrix of size \eqn{(256\times 256)}
#' @keywords datasets
#' @references Gonzalez, Rafael C. and Woods, Richard E. (2017) \emph{Digital Image Processing} (4th ed.). ISBN 0133356728.
#'
#' @source \href{http://sipi.usc.edu/database/?volume=misc}{USC SIPI Image Database}
#'
#' @examples
#' data(lena256)
#' image(lena256, col=gray((0:100)/100), axes=FALSE, main="lena256")
"lena256"
