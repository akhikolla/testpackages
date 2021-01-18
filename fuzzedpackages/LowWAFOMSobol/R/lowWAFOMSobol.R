##' Low WAFOM Sobol Sequence
##'
##' Description: R implementation of Low Walsh Figure of Merit (WAFOM)
##' Sequence based on Sobol Sequence.
##'
##' Porting to R by Mutsuo Saito.
##' The R version does not return coordinate value zero,
##' but returns value very near to zero, 2^-64.
##'
##'@section Acknowledgment:
##' The development of this code is partially supported
##' by JST CREST.
##'
##'@examples
##' srange <- lowWAFOMSobol.dimMinMax()
##' mrange <- lowWAFOMSobol.dimF2MinMax(srange[1])
##' points <- lowWAFOMSobol.points(dimR=srange[1], dimF2=mrange[1])
##' points <- lowWAFOMSobol.points(dimR=srange[1], dimF2=mrange[1], digitalShift=TRUE)
##'@section Reference:
##' * Shinsuke Mori,
##'   "Suuchi Sekibun no tameno QMC Ten Shuugou no Sekkei, Tansaku,
##'   oyobi sono Yuukousei",
##'   Master's Thesis, 2017,
##' * Ryuichi Ohori,
##'   "Efficient Quasi Monte Carlo Integration by Adjusting the
##'   Derivation-sensitivity Parameter of Walsh Figure of Merit",
##'   Master's Thesis, 2015.
##' * S. Harase and R. Ohori,
##'   "A search for extensible low-WAFOM point sets",
##'   arXiv preprint, arXiv:1309.7828, (2013),
##'   https://arxiv.org/abs/1309.7828.
##' * Harase, S. (2016).
##'   "A search for extensible low-WAFOM point sets",
##'   Monte Carlo Methods and Applications, 22(4), pp. 349-357, 2017.
##    doi:10.1515/mcma-2016-0119
##' * M. Matsumoto and R. Ohori,
##'   "Walsh Figure of Merit for Digital Nets: An Easy Measure
##'   for Higher Order Convergent QMC",
##'   Springer International Publishing, Cham, 2016, pp. 143-160.
##' * M. Matsumoto, M. Saito, and K. Matoba,
##'   "A computable figure of merit for quasi-Monte Carlo point sets",
##'   Mathematics of Computation, 83 (2014), pp. 1233-1250.
##' * S. Joe and F. Y. Kuo,
##'   "Constructing Sobol sequences with better two-dimensional projections",
##'   SIAM J. Sci. Comput. 30, 2635-2654 (2008).
##'
##'@name LowWAFOMSobol-package
##'@aliases LowWAFOMSobol-package LowWAFOMSobol
##'@docType package
##'@import Rcpp
##'@import RSQLite
##'@importFrom stats runif
##'@useDynLib LowWAFOMSobol
NULL

##' get minimum and maximum dimension number of Low WAFOM Niederreiter-Xing
##' Sequence
##'
##'@return supported minimum and maximum dimension number.
##'@export
lowWAFOMSobol.dimMinMax <- function() {
  fname <- "solw64.sqlite3"
  sql <- "select min(dimr) as min, max(dimr) as max from digitalnet;"
  drv <- dbDriver("SQLite")
  con <- dbConnect(drv, dbname = system.file("extdata",
                                             fname, package = "LowWAFOMSobol"))
  a <- dbGetQuery(con, sql)
  dbDisconnect(con)
  c(a[1, 1], a[1, 2])
}

##' get minimum and maximum F2 dimension number.
##'
##'@param dimR dimension.
##'@return supported minimum and maximum F2 dimension number
##'@export
lowWAFOMSobol.dimF2MinMax <- function(dimR) {
  fname <- "solw64.sqlite3"
  fmt <- "select min(dimf2) as min, max(dimf2) as max from digitalnet where dimr = %d;"
  sql <- sprintf(fmt, dimR)
  drv <- dbDriver("SQLite")
  con <- dbConnect(drv, dbname = system.file("extdata",
                                             fname, package = "LowWAFOMSobol"))
  a <- dbGetQuery(con, sql)
  dbDisconnect(con)
  c(a[1, 1], a[1, 2])
}

##' get points from Low WAFOM SobolSequence
##'
##' This R version does not returns coordinate value zero,
##' but returns value very near to zero, 2^-64.
##'@param dimR dimension.
##'@param dimF2 F2-dimension of each element.
##'@param digitalShift use digital shift or not.
##'@return matrix of points where every row contains dimR dimensional point.
##'@export
lowWAFOMSobol.points <- function(dimR,
                              dimF2 = 10,
                              digitalShift = FALSE) {
  smax = lowWAFOMSobol.dimMinMax()
  if (dimR < smax[1] || dimR > smax[2]) {
    stop(sprintf("dimR should be an integer %d <= dimR <= %d", smax[1], smax[2]))
  }
  mmax = lowWAFOMSobol.dimF2MinMax(dimR)
  if (dimF2 < mmax[1] || dimF2 > mmax[2]) {
    stop(sprintf("dimF2 should be an integer %d <= dimF2 <= %d", mmax[1], mmax[2]))
  }
  fname <- "solw64.sqlite3"
  fmt <- paste("select dimr, dimf2, wafom, tvalue, quote(data) as data from digitalnet ",
               "where dimr = %d and dimf2 = %d;")
  sql <- sprintf(fmt, dimR, dimF2)
  drv <- dbDriver("SQLite")
  con <- dbConnect(drv, dbname = system.file("extdata",
                                             fname, package = "LowWAFOMSobol"))
  df <- dbGetQuery(con, sql)
  dbDisconnect(con)
  if (digitalShift) {
    sv <- runif(2*dimR, min=-2^31, max=2^31-1)
  } else {
    sv <- numeric(1)
  }
  #  print(sv)
  count <- 2^dimF2
  return(rcppLowWAFOMSobolPoints(df, dimR, dimF2, count, sv))
}
