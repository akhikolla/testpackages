#' Wrapper to output the mutual information networks for subpopulations with size larger than a desired threshold.
#' 
#' @param dataMatrix data matrix with first column being the sample ID
#' @param subpopulationLabels the subpopulation labels
#' @param threshold the number of edges to draw for each subpopulation mutual information network

outputNetworks_topEdges_matrix<-function(dataMatrix, subpopulationLabels, threshold){
  
  SampleID<-dataMatrix[1,1]
  inputMatrix<-cbind(dataMatrix, subpopulationLabels)
  labels<-unique(subpopulationLabels)

  #plot and output networks
  for(i in 1:length(labels)){
    inputMatrix_subpopSubset<-subset(inputMatrix, inputMatrix[,ncol(inputMatrix)] == i)
    inputMatrix_subpopSubset_datapoints<-inputMatrix_subpopSubset[,c(2:(ncol(inputMatrix_subpopSubset)-1))]
    
    #plot networks
    png(paste(SampleID, '_Subpopulation', "0", i, '.png', sep=""))
    MINetworkPlot_topEdges(inputMatrix_subpopSubset_datapoints, threshold)
    title(main=paste("Subpopulation", i))
    dev.off()
    
    #output network matrix
    Networks<-MINetwork_matrix_topEdges(inputMatrix_subpopSubset_datapoints, threshold)
    
    save(Networks, file=paste0(SampleID, "_", i, "_subpop_Network.Rdata"))
  }
  
}