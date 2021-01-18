
#' Creates annotation matrix for the clades in aggregated format. The matrix contains average signals of each dimension for each clade in each sample
#' 
#' @param sampleIDs sampleID vector
#' @param topHubs number of top ranked genes to output for annotation; annotation is a concatenated list of top ranked genes.
#' @return Annotated clade matrix
#' @export

#annotation by network hub genes
annotateClades<-function(sampleIDs, topHubs){
  
  dataMatrix_filtered<-NULL
  newlabels<-NULL
  
  aggregateMatrix<-matrix(0,0,0)
  
  for(i in 1:length(sampleIDs)){
    sampleID<-sampleIDs[i]
    
    load(paste0(sampleID,"_dataMatrix_filtered.Rdata"))
    load(paste0(sampleID,"_new_subpopulations_Representative_subpops_filtered.Rdata"))
    
    inputMatrix<-as.data.frame(dataMatrix_filtered)
    data_agg<-aggregateData(inputMatrix,newlabels)
    
    aggregateMatrix<-rbind(aggregateMatrix,data_agg)
  }  
  
  
  aggregateRepNetworks<-matrix(0,0,0)
  mainDir <- getwd()
  
  for(i in 1:length(sampleIDs)){
    sampleID<-sampleIDs[i]
    subDir <- paste0("/", sampleID, "_CladeNetworks")
    
    setwd(paste0(mainDir, subDir))
    
    load(paste0(sampleID,"_Representative_Network.Rdata"))
    Networks<-as.data.frame(Networks)
    aggregateRepNetworks<-rbind(aggregateRepNetworks, Networks)
    setwd(mainDir)
  }
  
  subpopName<-matrix(0,0,1)
  for(i in 1:nrow(aggregateRepNetworks)){
    orderedGeneHubs<-sort(aggregateRepNetworks[i,], decreasing=TRUE)
    annotationGeneHubs<-orderedGeneHubs[1:topHubs]
    subpopName[i]<-paste(names(annotationGeneHubs), collapse="-")
    
  }
  
  annotationMatrix<-cbind(subpopName, rownames(aggregateRepNetworks))
  subpopulationRawNames<-paste0(aggregateMatrix[,2], "_",aggregateMatrix[,1])
  
  
  translateMatrix<-matrix(0,0,0)
  for(i in 1:length(subpopulationRawNames)){
    index<-which(subpopulationRawNames[i]==annotationMatrix[,2])
    translateMatrix[i]<-annotationMatrix[index,1]
  }
  
  aggregateMatrix_withAnnotation<-cbind(translateMatrix, aggregateMatrix)
  colnames(aggregateMatrix_withAnnotation)[1]<-"Annotation"
  
  write.table(aggregateMatrix_withAnnotation, file="aggregateMatrix_withAnnotation.xls", sep="\t", row.names=FALSE)
  return(aggregateMatrix_withAnnotation)
}
