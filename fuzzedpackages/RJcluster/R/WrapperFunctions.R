#' scaleRJ
#'
#' @param X Numeric matrix to be scaled
#' @param medians A logical, if TRUE the medians of the columns are used to scaled
#'
#' @return The scaled matrix
#' @export
#'
#' @examples
#' data = generateSimulationData()
#' Xscaled = scaleRJ(data$X)
scaleRJ = function( X, medians = FALSE )
{
  if ( !is.matrix( X ) )
  {
    stop( "Input must be in matrix form" )
  }

  X = scale_c( X, medians )

  return( X )
}
