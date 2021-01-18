#'The Max-Min Normalization
#'
#' This function normalizes the data using the max-min normalization
#' @param x the dataset.
#' @param margin data is normalized by row (margin = 1) or by column (margin = 2). The default is 2.
#' @author Wanwan Zheng
#'
#' @examples
#'
#' data(ozone)
#' scaled.data = normalization(ozone[,-1])
#' ozone.scale = data.frame(y = as.character(ozone[,1]), scaled.data[,-1])
#'
#'@name normalization
#'@export

normalization = function(x,margin = 2)
{
  as.data.frame(apply(x,margin,normalized))
}
