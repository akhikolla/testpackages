#' exdex: Estimation of the Extremal Index
#'
#' The extremal index \eqn{\theta} is a measure of the degree of local
#' dependence in the extremes of a stationary process.  The \emph{exdex}
#' package  performs frequentist inference about \eqn{\theta} using the
#' methodologies proposed in Northrop (2015), Berghaus and Bucher (2018),
#' Suveges (2007) and Suveges and Davison (2010).
#'
#' @details Functions to implement three estimators of the extremal index
#'   are provided, namely
#' \itemize{
#'   \item{\code{\link{spm}}: semiparametric maxima estimator, using block
#'     maxima: (Northrop, 2015; Berghaus and Bucher, 2018)}
#'   \item{\code{\link{kgaps}}: \eqn{K}-gaps estimator, using threshold
#'     interexceedance times (Suveges and Davison, 2010)}
#'   \item{\code{\link{iwls}}: iterated weighted least squares estimator,
#'     using threshold interexceedance times: (Suveges, 2007)}
#' }
#'
#'   See \code{vignette("exdex-vignette", package = "exdex")} for an
#'   overview of the package.
#' @references Berghaus, B., Bucher, A. (2018) Weak convergence of a pseudo
#' maximum likelihood estimator for the extremal index. \emph{Ann. Statist.}
#' \strong{46}(5), 2307-2335. \url{https://doi.org/10.1214/17-AOS1621}
#' @references Northrop, P. J. (2015) An efficient semiparametric maxima
#' estimator of the extremal index. \emph{Extremes} \strong{18}(4), 585-603.
#' \url{https://doi.org/10.1007/s10687-015-0221-5}
#' @references Suveges, M. (2007) Likelihood estimation of the extremal
#'   index. \emph{Extremes}, \strong{10}, 41-55.
#'   \url{https://doi.org/10.1007/s10687-007-0034-2}
#' @references Suveges, M. and Davison, A. C. (2010) Model
#'   misspecification in peaks over threshold analysis, \emph{The Annals of
#'   Applied Statistics}, \strong{4}(1), 203-221.
#'   \url{https://doi.org/10.1214/09-AOAS292}
#' @seealso \code{\link{spm}} for estimation of the extremal index
#'   \eqn{\theta} using a semiparametric maxima method.
#' @seealso \code{\link{kgaps}} for maximum likelihood estimation of the
#'   extremal index \eqn{\theta} using the \eqn{K}-gaps model.
#' @seealso \code{\link{iwls}}: iterated weighted least squares estimator.
#' @seealso \code{\link{newlyn}} and \code{\link{sp500}} for example datasets.
#' @docType package
#' @name exdex
#' @import methods
#' @importFrom graphics plot
#' @importFrom stats coef confint nobs vcov
#' @useDynLib exdex, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' Newlyn sea surges
#'
#' The vector \code{newlyn} contains 2894 maximum sea-surges measured at
#' Newlyn, Cornwall, UK over the period 1971-1976. The observations are
#' the maximum hourly sea-surge heights over contiguous 15-hour time
#' periods.
#' @format A vector of length 2894.
#' @source Coles, S.G. (1991) Modelling extreme multivariate events. PhD thesis,
#'   University of Sheffield, U.K.
#' @references Fawcett, L. and Walshaw, D. (2012) Estimating return levels from
#'   serially dependent extremes. \emph{Environmetrics}, \strong{23}(3),
#'   272-283.  \url{https://doi.org/10.1002/env.2133}
#' @references Northrop, P. J. (2015) An efficient semiparametric maxima
#'   estimator of the extremal index. \emph{Extremes}, \strong{18},
#'   585-603.  \url{https://doi.org/10.1007/s10687-015-0221-5}
"newlyn"

#' Daily log returns of the Standard and Poor (S&P) 500 index
#'
#' Daily log returns of the S&P 500 index, that is, the log of the ratio of
#' successive daily closing prices, from 3rd January 1990 to 9th October 2018.
#'
#' @format A vector of length 7250, created using \code{\link[zoo]{zoo}}
#'   with an "index" attribute giving the date of the corresponding negated
#'   log return.
#' @source Yahoo finance \url{https://finance.yahoo.com/quote/^SPX/history/}
"sp500"
