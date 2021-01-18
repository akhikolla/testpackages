##' Finite Mixture of Hidden Markov Models for accelerometer data
##'
##' \tabular{ll}{
##'   Package: \tab MHMM\cr
##'   Type: \tab Package\cr
##'   Version: \tab 1.0.0\cr
##'   Date: \tab 2020-03-20\cr
##'   License: \tab GPL-2\cr
##'   LazyLoad: \tab yes\cr
##' }
##'
##'
##' @name MHMM-package
##' @aliases MHMM
##' @rdname MHMM-package
##' @docType package
##' @keywords package
##' @import Rcpp
##' @import parallel
##' @import methods
##' @import ggplot2
##' @import reshape2
##' @import gridExtra
##' @importFrom stats density dgamma na.omit quantile rgamma rnorm rpois runif var
##' @importFrom grDevices hcl
##' @useDynLib MHMM
##' @references  Du Roy de Chaumaray, M. and Marbac, M. and Navarro, F. (2019). Mixture of hidden Markov models for accelerometer data. arXiv preprint arXiv:1906.01547
##' @examples
##' data(accelero)
##' # To make the estimation <5
##' res <- mhmm(accelero, K = 2, M = 4, nbcores = 1, nbinit = 5, iterSmall = 2)
##' plot(res, 1)
##' 
##'  \donttest{
##' data(accelero)
##' # It is better to increase the number of random initializations
##' res <- mhmm(accelero, K = 2, M = 4, nbcores = 1)
##' plot(res, 1)
##' }
NULL


##' Accelerometer data
##'
##' Accelerometer data measured each 5 minutes on three subjects
##'
##'
##'
##' @references  Huang, Q., Cohen, D., Komarzynski, S., Li, X.-M., Innominato, P., Lévi, F., and Finkenstädt, B. (2018b). Hidden markov models for monitoring circadian rhythmicity in telemetric activity data. Journal of The Royal Society Interface, 15(139):20170885
##' @name accelero
##' @docType data
##' @keywords datasets
##'
##' @examples
##'   data(accelero)
NULL
