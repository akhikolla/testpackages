##' Sobol Sequence
##'
##' R implementation of S. Joe and F. Y. Kuo,
##' "Constructing Sobol sequences with better two-dimensional projections",
##' SIAM J. Sci. Comput. 30, 2635-2654 (2008).
##'
##' The implementation is based on the data file new-joe-kuo-6.21201
##' <http://web.maths.unsw.edu.au/~fkuo/sobol/>.
##'
##' Porting to R by Mutsuo Saito.
##' The R version does not returns cordinate value zero,
##' but returns value very near to zero, 2^-64.
##'
##'@section Acknowledgments:
##' I, Mutsuo Saito, wish to thank Frances Kuo and Stephen Joe for their research,
##' and agreement to use thier source code.
##'
##' The development of this R code is partially supported by JST CREST.
##'
##'@examples
##' srange <- sobolSequence.dimMinMax()
##' mrange <- sobolSequence.dimF2MinMax(srange[1])
##' points <- sobolSequence.points(dimR=srange[1], dimF2=mrange[1], count=10000)
##' points <- sobolSequence.points(dimR=srange[1], dimF2=mrange[1], count=10000,
##'                                digitalShift=TRUE)
##'@section Reference:
##' S. Joe and F. Y. Kuo,
##' "Constructing Sobol sequences with better two-dimensional projections",
##' SIAM J. Sci. Comput. 30, 2635-2654 (2008).
##'
##'@name SobolSequence-package
##'@aliases SobolSequence-package sobolsequence
##'@docType package
##'@import Rcpp
##'@importFrom stats runif
##'@useDynLib SobolSequence
NULL

##' get minimum and maximum dimension number of Sobol Sequence
##'
##'@return supportd minimum and maximum dimension number.
##'@export
sobolSequence.dimMinMax <- function() {
    return(c(2, 21201))
}

##' get minimum and maximum F2 dimension number.
##'
##'@param dimR dimention.
##'@return supportd minimum and maximum F2 dimension number.
##'@export
sobolSequence.dimF2MinMax <- function(dimR) {
    return(c(10, 31))
}

##' get points from SobolSequence
##'
##' This R version does not returns cordinate value zero,
##' but returns value very near to zero, 2^-64.
##'@param dimR dimention.
##'@param dimF2 F2-dimention of each element.
##'@param count number of points.
##'@param digitalShift use digital shift or not.
##'@return matrix of points where every row contains dimR dimensional point.
##'@export
sobolSequence.points <- function(dimR,
                              dimF2 = 10,
                              count,
                              digitalShift = FALSE) {
  smax = sobolSequence.dimMinMax()
  if (dimR < smax[1] || dimR > smax[2]) {
    stop(sprintf("dimR should be an integer %d <= dimR <= %d", smax[1], smax[2]))
  }
  if (missing(dimF2)) {
    dimF2 = max(dimF2, ceiling(log2(count)))
  }
  mmax = sobolSequence.dimF2MinMax(dimR)
  if (dimF2 < mmax[1] || dimF2 > mmax[2]) {
    stop(sprintf("dimF2 should be an integer %d <= dimF2 <= %d", mmax[1], mmax[2]))
  }
  fname = system.file("extdata",
                       "new-joe-kuo-6.21201",
                       package = "SobolSequence")
  if (digitalShift) {
    sv <- runif(2*dimR, min=-2^31, max=2^31-1)
  } else {
    sv <- numeric(1)
  }
#  print(sv)
  return(rcppSobolPoints(fname, dimR, dimF2, count, sv))
}
