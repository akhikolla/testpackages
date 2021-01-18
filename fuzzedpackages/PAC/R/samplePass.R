#' Run PAC for Specified Samples
#' 
#' A wrapper to run PAC and output subpopulation mutual information networks. Please use the PAC function itself for individual samples or if the MAN step is not needed.
#' 
#' 
#' @param sampleIDs sampleID vector
#' @param dim_subset a string vector of string names to subset the data columns for PAC; set to NULL to use all columns
#' @param hyperrectangles number of hyperrectangles to learn for each sample
#' @param num_PACSupop number of subpopulations to output for each sample using PAC
#' @param max.iter postprocessing kmeans iterations
#' @param num_networkEdge a threshold on the number of edges to output for each subpopulation mutual information network
#' @export
#' @examples 
#' Please see vignette for the PAC-MAN tutorial

samplePass<-function(sampleIDs, dim_subset, hyperrectangles,num_PACSupop, max.iter,num_networkEdge){
  for(i in 1:length(sampleIDs)){
    
    load(paste0(sampleIDs[i], ".Rdata"))
    inputMatrix<-get(sampleIDs[i])
    
    #subset for dimensions involved in clustering of subpopulations
    if(!is.null(dim_subset)){
      inputMatrix_dimsubset<-inputMatrix[,dim_subset]
      
      subpopulationLabels = PAC(inputMatrix_dimsubset, num_PACSupop, maxlevel = hyperrectangles, method = "dsp", max.iter = max.iter)
      save(subpopulationLabels, file=paste0(sampleIDs[i], "_subpopulationLabels.Rdata"))
      
      inputMatrix_withSampleName<-cbind(sampleIDs[i], as.data.frame(inputMatrix))
      data_agg<-aggregateData(inputMatrix_withSampleName,subpopulationLabels)
      save(data_agg, file=paste0(sampleIDs[i], "_data_agg.Rdata"))
      save(inputMatrix_withSampleName, file=paste0(sampleIDs[i], "_dataMatrix.Rdata"))
      
      mainDir <- getwd()
      subDir <- paste0("/", sampleIDs[i], "_SampleSubpopulationMatrixNetworks")
      dir.create(file.path(mainDir, subDir))
      setwd(paste0(mainDir, subDir))
      inputMatrix_dimsubset_withSampleName<-cbind(sampleIDs[i], as.data.frame(inputMatrix_dimsubset))
      
      outputNetworks_topEdges_matrix(inputMatrix_dimsubset_withSampleName, subpopulationLabels, threshold=(2*num_networkEdge))
      setwd(mainDir)
      
    }else{
    
    subpopulationLabels = PAC(inputMatrix, num_PACSupop, maxlevel = hyperrectangles, method = "dsp", max.iter = max.iter)
    save(subpopulationLabels, file=paste0(sampleIDs[i], "_subpopulationLabels.Rdata"))
    
    inputMatrix_withSampleName<-cbind(sampleIDs[i], as.data.frame(inputMatrix))
    data_agg<-aggregateData(inputMatrix_withSampleName,subpopulationLabels)
    save(data_agg, file=paste0(sampleIDs[i], "_data_agg.Rdata"))
    save(inputMatrix_withSampleName, file=paste0(sampleIDs[i], "_dataMatrix.Rdata"))
    
    mainDir <- getwd()
    subDir <- paste0("/", sampleIDs[i], "_SampleSubpopulationMatrixNetworks")
    dir.create(file.path(mainDir, subDir))
    setwd(paste0(mainDir, subDir))
    outputNetworks_topEdges_matrix(inputMatrix_withSampleName, subpopulationLabels, threshold=(2*num_networkEdge))
    setwd(mainDir)
    
    }
  }
}