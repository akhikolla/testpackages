#' Simulated test data.
#' 
#' @description A dataset containing time series for 100 variates with 1000 data points for each variate.
#' The dataset contains 5 changepoints with each one being shared independently by 10 variates.
#' The observations prior to each changepoint are IID Gaussian, distributed with unit variance and random mean drawn from N(0,1) Gaussians.
#' The mean after each changepoint, for each variable that is selected to have a change in mean, is changed by size=1 at each changepoint location with the sign of change chosen uniformly at random.
#' 
#' @name testdata
#' 
#' @docType data
#' 
#' @keywords dataset
#' 
#' @usage data(testdata)
#' 
#' @rdname testdata-data
#' 
#' @format A matrix with 100 rows and 2000 columns.
#' 
#' @references Hahn, G., Fearnhead, P., Eckley, I.A. (2020). Fast computation of a projection direction for multivariate changepoint detection. Stat Comput.
#' 
NULL
