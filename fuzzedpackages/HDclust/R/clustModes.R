#' @title Hierarchical clustering of density modes
#'
#' @description This function performs hierarchical clustering of density modes found
#' by \code{hmmvbFindModes()}. 
#' @param dist.args A list with arguments to \code{dist} method from package \code{stats}. It is used to 
#' specify a distance metric for the density modes. Euclidian distance is used by default. See \code{dist} for details.
#' @param hclust.args A list with arguments to \code{hclust} method from package \code{stats}. It is used to 
#' specify a linkage method for hierarchical clustering of density modes. Complete linkage is used by default.
#' See \code{hclust} for details.
#' @param cutree.args A list with arguments to \code{cutree} method from package \code{stats}.
#' @param modes An object of class 'HMMVBclust' returned by \code{hmmvbFindModes()}. 
#' @return An object of class 'HMMVBclust' with new cluster labels and cluster sizes. Note that coordinates of modes after merging are not calculated and 
#' \code{clustParam} field is empty.
#' @seealso \code{\link{hmmvbClust}}, \code{\link{hmmvbFindModes}}
#' @examples
#' Vb <- vb(1, dim=4, numst=2)
#' set.seed(12345)
#' hmmvb <- hmmvbTrain(unique(iris[,1:4]), VbStructure=Vb)
#' modes <- hmmvbFindModes(unique(iris[,1:4]), model=hmmvb)
#'
#' # default mode clustering
#' merged <- clustModes(modes, cutree.args=list(h=1.0))
#' 
#' # mode clustering using Manhattan distance
#' merged <- clustModes(modes, dist.args=list(method="manhattan"), cutree.args=list(h=1.0))
#' 
#' # mode clustering using single linkage
#' merged <- clustModes(modes, hclust.args=list(method="single"), cutree.args=list(h=1.0))
clustModes <- function(modes, cutree.args, hclust.args=NULL, dist.args=NULL){
  
  clustParam <- getClustParam(modes)
  
  normModes <- t(t(clustParam$mode) / clustParam$sigma)
  
  if (is.null(dist.args))
    dist.args = list(normModes)
  else
    dist.args = c(list(normModes), dist.args)
  
  d <- do.call(dist, dist.args)
  
  if (is.null(hclust.args))
    hclust.args = list(d)
  else
    hclust.args = c(list(d), hclust.args)
  
  
  hc <- do.call(hclust, hclust.args)
  
  cutree.args = c(list(hc), cutree.args)
  
  newModeid <- do.call(cutree, cutree.args)
  #print(newModeid)
  newClsid <- newModeid[as.factor(getClsid(modes))]
  newSize <- rep(0, max(newModeid))
  
  for (i in 1:max(newModeid))
    newSize[i] = length(which(newClsid == i))
  
  summary(hc)
  plot(as.dendrogram(hc), ylab="Height")
  
  return(new("HMMVBclust", data=modes@data, clsid=newClsid, size=newSize))
  #return(list(data=modes@data, clsId=newClsId, size=newSize))
  
}