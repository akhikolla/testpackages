#' Calculates the within cluster spread 
#' @param sampleIDs A vector of sample names.
#' @param dim_subset a string vector of string names to subset the data columns for PAC; set to NULL to use all columns.
#' @param SubpopSizeFilter threshold to filter out very small clusters with too few points; these very small subpopulations may not be outliers and not biologically relevant.
#' 
#' @return Returns the sample within cluster spread
#' 

recordWithinClusterSpread<-function(sampleIDs,dim_subset=NULL,SubpopSizeFilter){
  sample_mean_dist<-vector()
  
  for(j in 1:length(sampleIDs)){
    newer_subpopulationLabels <- NULL
    inputMatrix_withSampleName <- NULL
    
    sampleID <- sampleIDs[j]
    load(paste0(sampleID, "_new_subpopulations_Representative_subpops.Rdata"))
    load(paste0(sampleID, "_dataMatrix.Rdata"))
    
    newlabels <- newer_subpopulationLabels
    filteredClades <- names(table(newlabels))[table(newlabels) > SubpopSizeFilter]
    whetherToKeep <- newlabels %in% filteredClades
    newlabels <- newlabels[whetherToKeep]
    
    dataMatrix_filtered <- inputMatrix_withSampleName[whetherToKeep,]
    
    if (!is.null(dim_subset)) {
      dataMatrix_filtered <- cbind(dataMatrix_filtered[,1], dataMatrix_filtered[, dim_subset])
    }
    
    dist_record<-vector()
    
    subpop_proportions<-table(newlabels)/length(newlabels)
    
    for(i in 1:length(unique(newlabels))){
      
      newlabels_iter<-unique(newlabels)[i]
      
      newer_subpopulationLabels_index<-which(newlabels == newlabels_iter)
      
      sample_cladeSubset<-dataMatrix_filtered[newer_subpopulationLabels_index,-1]
      
      matrix_centroid<-colSums(sample_cladeSubset)/nrow(sample_cladeSubset)
      
      #weighted average spread by subpopulation proportion
      
      distanceFromCentroid<-sweep(sample_cladeSubset,2, matrix_centroid)
      dist_record[i]<-mean(sqrt(rowSums( (distanceFromCentroid)^2 ) ))*(subpop_proportions[newlabels_iter])
      
    }
    
    sample_mean_dist[j]<-mean(dist_record)
    
  }
  return(sample_mean_dist)
}
