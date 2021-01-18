#' @title Parameters for MBM clustering algorithm.
#' @name clustControl
#' @description This function creates a list with parameters for Modal Baum-Welch (MBW)
#' clustering algorithm used as an argument for \code{hmmvbClust}.
#' @param minSize Minimum cluster size. Clusters that contain the number of data points
#' smaller than \code{minSize} are merged to the closest big cluster.
#' @param modeTh Distance parameter that controls mode merging. Larger values 
#' promote merging of different clusters. 
#' @param useL1norm A logical value indicating whether or not L1 norm will be 
#' used to calculate the distance.
#' @param getlikelh A logical value indicating whether or not to calculate the
#' loglikelihood for every data point.
#' @return The named list with parameters.
#' @seealso \code{\link{hmmvbTrain}}
#' @examples 
#' # avoid clusters of size < 60
#' Vb <- vb(1, dim=4, numst=2)
#' set.seed(12345)
#' hmmvb <- hmmvbTrain(iris[,1:4], VbStructure=Vb)
#' clust <- hmmvbClust(iris[,1:4], model=hmmvb, control=clustControl(minSize=60))
#' show(clust)
clustControl <- function(minSize=1, modeTh=0.01, useL1norm=FALSE, getlikelh=FALSE){
  if ((!is.numeric(minSize)) || (length(minSize)!=1) || (minSize <= 0) || (minSize != as.integer(minSize)))
    stop('Invalid minSize argument provided! minSize should be a scalar positive integer\n')
  
  if ((!is.numeric(modeTh)) || (length(modeTh)!=1))
    stop('Invalid modeTh argument provided! modeTh should be a scalar double\n')
  
  if ((!is.logical(useL1norm)) || (length(useL1norm)!=1))
    stop('Invalid useL1norm argument provided! useL1norm should be a logical\n')
    
  if ((!is.logical(getlikelh)) || (length(getlikelh)!=1))
  stop('Invalid getlikelh argument provided! getlikelh should be a logical\n')
  
  return(list(minSize=minSize, modeTh=modeTh, useL1norm=useL1norm, getlikelh=getlikelh))
}
