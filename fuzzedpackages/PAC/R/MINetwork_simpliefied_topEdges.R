#' Outputs the vectorized summary of a network based on the number of edges connected to a node
#' @import parmigene
#' @importFrom infotheo discretize mutinformation
#' @param dataMatrix data matrix
#' @param threshold the number of edges to draw for each subpopulation mutual information network

MINetwork_simplified_topEdges<-function(dataMatrix, threshold){
  
  dat<-discretize(dataMatrix)
  dataMatrix_mi<-mutinformation(dat, method="emp")
  
  dataMatrix_net<-mrnet(dataMatrix_mi)
  
  filterLevel<-dataMatrix_net[order(dataMatrix_net, decreasing=TRUE)][threshold]
  dataMatrix_net[dataMatrix_net<=filterLevel]<-0
  dataMatrix_net[dataMatrix_net>=filterLevel]<-1
  return(colSums(dataMatrix_net))
  
}