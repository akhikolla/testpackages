if ( !isGeneric("meanDifference") ) {
  setGeneric("meanDifference", function(x, ...)
    standardGeneric("meanDifference"))
}
#' Calculate mean difference between two datasets
#'
#' @description
#' Calculate the mean difference between two datasets as suggested by 
#' Wang \emph{et al.} (2012). 
#'
#' @param x,y Objects of class \code{RasterLayer} or \code{numeric}. 
#'
#' @return
#' The mean difference between the two inputs either as \code{RasterLayer} or 
#' \code{numeric}.
#'  
#' @source 
#' Wang \emph{et al.} (2012) Impact of sensor degradation on the MODIS NDVI time 
#' series. Remote Sensing of Environment 119, 55-61, 
#' doi:\href{https://doi.org/10.1016/j.rse.2011.12.001}{10.1016/j.rse.2011.12.001}.
#' 
#' Detsch \emph{et al.} (2016) A Comparative Study of Cross-Product NDVI 
#' Dynamics in the Kilimanjaro Region - A Matter of Sensor, Degradation 
#' Calibration, and Significance. Remote Sensing 8(2), 159, 
#' doi:\href{http://dx.doi.org/10.3390/rs8020159}{10.3390/rs8020159}.
#'     
#' @examples
#' x <- 1:10
#' y <- 2:11
#' meanDifference(x, y)
#'
#' @export meanDifference
#' @name meanDifference

################################################################################
### function using 'RasterLayer' objects ---------------------------------------
#' @aliases meanDifference,RasterLayer-method
#' @rdname meanDifference
setMethod("meanDifference",
          signature(x = 'RasterLayer'),
          function(x, y) {
  
  ## merge values
  num_x <- x[]
  num_y <- y[]
  
  if (length(num_x) != length(num_y))
    cat("Warning: elements 'x' and 'y' have unequal number of cells.\n")
  
  num_xy <- cbind(num_x, num_y)
  num_xy <- num_xy[stats::complete.cases(num_xy), ]
  
  ## calculate mean difference
  mean(num_xy[, 1] - num_xy[, 2])
})

################################################################################
### function using 'numeric' vectors -------------------------------------------
#' @aliases meanDifference,numeric-method
#' @rdname meanDifference
setMethod("meanDifference",
          signature(x = "numeric"),
          function(x, y) {

  if (length(x) != length(y))
    cat("Warning: elements 'x' and 'y' are of unequal length.\n")
  
  ## merge values
  num_xy <- cbind(x, y)
  num_xy <- num_xy[stats::complete.cases(num_xy), ]
  
  ## calculate mean difference
  mean(num_xy[, 1] - num_xy[, 2])
})
