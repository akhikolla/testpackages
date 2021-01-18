#' @import parmigene
#' @importFrom infotheo discretize mutinformation
#' @importFrom igraph graph.adjacency simplify igraph_opt layout.circle
NULL

#' Plots mutual information network (mrnet algorithm) connection using the parmigene package. Mutual information calculated with infotheo package.
#' @param dataMatrix data matrix
#' @param threshold the maximum number of edges to draw for each subpopulation mutual information network
#'
#' @export
#' 

MINetworkPlot_topEdges<-function(dataMatrix, threshold){

  dat<-discretize(dataMatrix)
  dataMatrix_mi<-mutinformation(dat, method="emp")
  
  dataMatrix_net<-mrnet(dataMatrix_mi)
  
  filterLevel<-dataMatrix_net[order(dataMatrix_net, decreasing=TRUE)][threshold]
  dataMatrix_net[dataMatrix_net<filterLevel]<-0
  
  g_dataMatrix_net <- graph.adjacency(dataMatrix_net, weighted=TRUE, mode="undirected")
  g_dataMatrix_net<-simplify(g_dataMatrix_net, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb=igraph_opt("edge.attr.comb"))
  networkPlot<-plot(g_dataMatrix_net, vertex.color=rainbow(ncol(dataMatrix)), edge.width=3,layout=layout.circle(g_dataMatrix_net))
}
