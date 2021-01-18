

#' Simulated data.
#'
#' @description A dataset containing time series for 100 variates with 500 data points
#' for each variate. The data contains 5 most recent changepoints (mrc's)
#' with each mrc being shared independently by 20 variates. The observations prior to
#' the mrcs are IID Gaussian distributed with unit variance and mean
#' \eqn{\mu \sim N(0,2^{2})}. The mean
#' after the mrc is \eqn{\mu \pm 0.5} with the sign of the change
#' chosen uniformly at random for each variate. For more details regarding
#' simulated testing of \link{mrc} see Bardwell, Eckley, Fearnhead and
#' Smith 2016, pages 5--7.
#' 
#' @name mrcexample
#' 
#' @docType data
#'
#' @keywords datasets
#'
#' @usage data(mrcexample)
#'
#' @rdname mrcexample-data
#'
#' @format A matrix with 100 rows and 500 columns.
#'
#' @references \insertRef{doi:10.1080/00401706.2018.1438926}{changepoint.mv}
#' 
NULL




