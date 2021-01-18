#' @import parmigene
#' @importFrom infotheo discretize mutinformation
NULL

#' Mutual information network connection matrix generation (mrnet algorithm) using the parmigene package. Mutual information calculated with infotheo package.
#' 
#' @param dataMatrix data matrix
#' @param threshold the number of edges to draw for each subpopulation mutual information network
#' @return the mutual information network connection matrix with top edges


MINetwork_matrix_topEdges<-function(dataMatrix, threshold){
  
  dat<-discretize(dataMatrix)
  dataMatrix_mi<-mutinformation(dat, method="emp")
  
  dataMatrix_net<-mrnet(dataMatrix_mi)
  
  filterLevel<-dataMatrix_net[order(dataMatrix_net, decreasing=TRUE)][threshold]
  dataMatrix_net[dataMatrix_net<=filterLevel]<-0
  dataMatrix_net[dataMatrix_net>=filterLevel]<-1
  return(dataMatrix_net)
  
}