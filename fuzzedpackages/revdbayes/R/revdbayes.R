#' revdbayes: Ratio-of-Uniforms Sampling for Bayesian Extreme Value Analysis
#'
#' Uses the multivariate generalized ratio-of-uniforms method to simulate
#' random samples from the posterior distributions commonly encountered in
#' Bayesian extreme value analyses.
#'
#' @details The main functions in the revdbayes package are \code{\link{rpost}}
#'   and \code{\link{rpost_rcpp}}, which simulate random samples from the
#'   posterior distribution of extreme value model parameters using the
#'   functions \code{\link[rust]{ru}} and \code{\link[rust]{ru_rcpp}}
#'   from the rust package, respectively. The user chooses the extreme value
#'   model, the prior density for the parameters and provides the data.
#'   There are options to improve the probability of acceptance of the
#'   ratio-of-uniforms algorithm by working with transformation of the model
#'   parameters.
#'
#'   The function \code{\link{kgaps_post}} simulates from the posterior
#'   distribution of the extremal index \eqn{\theta} based on the
#'   K-gaps model for threshold interexceedance times of Suveges and Davison
#'   (2010).  See also Attalides (2015).
#'
#'   See \code{vignette("revdbayes-a-vignette", package = "revdbayes")} for an
#'   overview of the package and
#'   \code{vignette("revdbayes-b-using-rcpp-vignette", package = "revdbayes")}
#'    for an illustration of the improvements in efficiency produced using
#'    the Rcpp package.
#'    See `vignette("revdbayes-c-predictive-vignette", package = "revdbayes")`
#'    for an outline of how to use revdbayes to perform posterior predictive
#'    extreme value inference.
#' @references Northrop, P. J. (2016). rust: Ratio-of-Uniforms Simulation with
#'   Transformation. R package version 1.2.2.
#'   \url{https://cran.r-project.org/package=rust}.
#' @references Suveges, M. and Davison, A. C. (2010) Model
#'   misspecification in peaks over threshold analysis, \emph{The Annals of
#'   Applied Statistics}, \strong{4}(1), 203-221.
#'   \url{https://doi.org/10.1214/09-AOAS292}
#' @references Attalides, N. (2015) Threshold-based extreme value modelling,
#'   PhD thesis, University College London.
#' @seealso \code{\link{set_prior}} to set a prior density for extreme value
#'   parameters.
#' @seealso \code{\link{rpost}} and \code{\link{rpost_rcpp}} to perform
#'   ratio-of-uniforms sampling from an extreme value posterior distribution.
#' @seealso \code{\link{kgaps_post}} to sample from the posterior distribution
#'   for the extremal index.
#' @seealso The \code{\link[rust]{ru}} and \code{\link[rust]{ru_rcpp}}
#'   functions in the \code{\link{rust}} package for details of the arguments
#'   that can be passed to \code{ru} via \code{rpost} and for the form of the
#'   object (of class "evpost") returned from \code{rpost}, which has the same
#'   structure as an object (of class "ru") returned by \code{ru} and
#'   \code{ru_rcpp}.
#' @docType package
#' @name revdbayes
#' @importFrom graphics plot
#' @importFrom stats predict
#' @importFrom bayesplot ppc_stat
#' @importFrom bayesplot ppc_stat_2d
#' @importFrom bayesplot ppc_dens_overlay
#' @importFrom bayesplot ppc_ecdf_overlay
#' @importFrom bayesplot ppc_dens
#' @importFrom bayesplot ppc_hist
#' @importFrom bayesplot ppc_boxplot
#' @importFrom bayesplot ppc_intervals
#' @importFrom bayesplot ppc_ribbon
#' @importFrom bayesplot pp_check
#' @importFrom bayesplot mcmc_areas
#' @importFrom bayesplot mcmc_intervals
#' @importFrom bayesplot mcmc_dens
#' @importFrom bayesplot mcmc_hist
#' @useDynLib revdbayes, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' Annual Maximum Sea Levels at Port Pirie, South Australia
#'
#' A numeric vector of length 65 containing annual maximum sea levels,
#' in metres, from 1923 to 1987 at Port Pirie, South Australia.
#'
#' @format A numeric vector containing 65 observations.
#' @source Coles, S. G. (2001) \emph{An Introduction to Statistical Modelling
#'   of Extreme Values}. London: Springer.
#'   \url{https://doi.org/10.1007/978-1-4471-3675-0}
"portpirie"

#' Annual Maximum Temperatures at Oxford
#'
#' A numeric vector containing annual maximum temperatures, in degrees
#' Fahrenheit, from 1901 to 1980 at Oxford, England.
#'
#'@format A vector containing 80 observations.
#'@source Tabony, R. C. (1983) Extreme value analysis in meteorology.
#'  \emph{The Meteorological Magazine}, \strong{112}, 77-98.
"oxford"

#' Daily Aggregate Rainfall
#'
#' A numeric vector of length 20820 containing daily aggregate rainfall
#' observations, in millimetres, recorded at a rain gauge in England
#' over a 57 year period, beginning on a leap year. Three of these years
#' contain only missing values.
#'
#' @format A vector containing 20820 observations.
#' @source Unknown
"rainfall"

#' Storm peak significant wave heights from the Gulf of Mexico
#'
#' A numeric vector containing 315 hindcasts of storm peak significant wave
#' heights, metres, from 1900 to 2005 at an unnamed location in the Gulf
#' of Mexico.
#'
#'@format A vector containing 315 observations.
#'@source Oceanweather Inc. (2005) GOMOS -- Gulf of Mexico hindcast study.
#' @references Northrop, P. J., Attalides, N. and Jonathan, P. (2017)
#'   Cross-validatory extreme value threshold selection and uncertainty
#'   with application to ocean storm severity.
#'   \emph{Journal of the Royal Statistical Society Series C: Applied
#'   Statistics}, \strong{66}(1), 93-120.
#'   \url{https://doi.org/10.1111/rssc.12159}
"gom"

#' Largest Sea Levels in Venice
#'
#' The \code{venice} data frame has 51 rows and 10 columns. The jth column
#' contains the jth largest sea levels in Venice, for the years 1931-1981. Only
#' the largest six measurements are available for the year 1935; the
#' corresponding row contains four missing values. The years for each set of
#' measurements are given as row names.
#'
#' @format A data frame with 51 rows and 10 columns.
#' @source Smith, R. L. (1986) Extreme value theory based on the \emph{r}
#'   largest annual events. \emph{Journal of Hydrology}, \strong{86}, 27-43.
#'   \url{https://doi.org/10.1016/0022-1694(86)90004-1}
#' @references Coles, S. G. (2001) \emph{An Introduction to Statistical
#'   Modelling of Extreme Values}. London: Springer.
#'   \url{https://doi.org/10.1007/978-1-4471-3675-0}
"venice"

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
