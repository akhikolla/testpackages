#' Outputs the representative/clade networks (plots and summary vectors) for subpopulations with size larger than a desired threshold. Saves the networks and the data matrices without the smaller subpopulations.
#' 
#' @param dataMatrix data matrix with first column being the sample ID
#' @param subpopulationLabels the subpopulation labels
#' @param threshold the number of edges to draw for each subpopulation mutual information network


outputRepresentativeNetworks_topEdges<-function(dataMatrix, subpopulationLabels, threshold){
  
  SampleID<-dataMatrix[1,1]
  inputMatrix<-cbind(dataMatrix, subpopulationLabels)
  labels<-unique(subpopulationLabels)
  
  #plot representative networks
  for(i in 1:length(labels)){
    inputMatrix_subpopSubset<-subset(inputMatrix, inputMatrix[,ncol(inputMatrix)] == labels[i])
    inputMatrix_subpopSubset_datapoints<-inputMatrix_subpopSubset[,c(2:(ncol(inputMatrix_subpopSubset)-1))]
    png(paste(SampleID, '_Representative_Subpopulation', "0", i, '.png', sep=""))
    MINetworkPlot_topEdges(inputMatrix_subpopSubset_datapoints, threshold)
    title(main=paste("Subpopulation", labels[i]))
    dev.off()
  }
  
  #output summarized network structures based on edge counts
  Networks<-matrix(0, length(labels), ncol(inputMatrix)-2)
  colnames(Networks)<-colnames(inputMatrix[,2:(ncol(inputMatrix)-1)])
  rownames(Networks)<-paste(SampleID, "_", labels, sep="")
  
  for(i in 1:length(labels)){
    inputMatrix_subpopSubset<-subset(inputMatrix, inputMatrix[,ncol(inputMatrix)] == labels[i])
    inputMatrix_subpopSubset_datapoints<-inputMatrix_subpopSubset[,c(2:(ncol(inputMatrix_subpopSubset)-1))]
    Networks[i,]<-MINetwork_simplified_topEdges(inputMatrix_subpopSubset_datapoints, threshold)
  }
  
  save(Networks, file=paste0(SampleID, "_Representative_Network.Rdata"))
}
