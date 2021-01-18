#' IMaGES Test Data
#'
#' Test data to be used with IMaGES. Included is
#' a list called \code{IMData}, a list containing three
#' data files that can be passed into IMaGES.
#' 
#' @name IMData
#'
#' @docType data
#'
#' @usage data(IMData)
#'
#' @format Objects of class \code{"list"}.
#'
#' @keywords data
#'
#'
#' @examples
#' data(IMData)
#' 
#' result <- IMaGES(matrices=IMData)
#' 
#' plotAll(result)
#' 
#' 
#' 
"IMData"

#' IMaGES Test Data
#'
#' Test data to be used with IMaGES. Included is
#' a list called \code{IMTrue}, a list containing the original
#' graph structures for each dataset in \code{IMTrue}.
#' 
#' @name IMTrue
#'
#' @docType data
#'
#' @usage data(IMData)
#'
#' @format Objects of class \code{"list"}.
#'
#' @keywords data
#'
#'
#' @examples
#' data(IMData)
#' data(IMTrue)
#' 
#' result <- IMaGES(matrices=IMData)
#' 
#' plotAll(result)
#' 
#' par(mfrow=c(1,length(IMTrue)))
#' 
#' for (i in 1:length(IMTrue)) {
#'   graph::plot(IMTrue[[i]])
#' }
#' 
#' 
"IMTrue"