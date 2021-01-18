
#' Runs elbow point analysis to find the practical optimal number of clades to output. Outputs the average within sample cluster spread for all samples and the elbow point analysis plot with loess line fitted through the results.
#' @param ks Vector that is a sequence of clade sizes.
#' @param sampleIDs A vector of sample names.
#' @param dim_subset a string vector of string names to subset the data columns for PAC; set to NULL to use all columns.
#' @param num_PACSupop Number of PAC subpopulation explored in each sample.
#' @param smallSubpopCutoff Cutoff of minor subpopulation not used in multiple alignments of networks
#' @param expressionGroupClamp clamps the subpopulations into desired number of expression groups for assigning small subpopulations into larger groups or their own groups.
#' @param SubpopSizeFilter threshold to filter out very small clusters with too few points in the calculation of cluster spreads; these very small subpopulations may be outliers and not biologically relevant.
#'
#' @export
#' 

runElbowPointAnalysis<-function(ks, sampleIDs, dim_subset, num_PACSupop, smallSubpopCutoff, expressionGroupClamp, SubpopSizeFilter){

  sample_mean_dist_MatrixRecord<-matrix(0, nrow=length(sampleIDs), ncol=length(ks))
  
  for(k_index in 1:length(ks)){
    
    print(k_index)
    clades_network_only<-MAN(sampleIDs, num_PACSupop=num_PACSupop, smallSubpopCutoff=smallSubpopCutoff, k_clades=ks[k_index])
    refineSubpopulationLabels(sampleIDs,dim_subset=dim_subset, clades_network_only, expressionGroupClamp=expressionGroupClamp)
    
    sample_mean_dist_MatrixRecord_update<-recordWithinClusterSpread(sampleIDs, dim_subset = dim_subset, SubpopSizeFilter = SubpopSizeFilter)
    sample_mean_dist_MatrixRecord[,k_index]<- sample_mean_dist_MatrixRecord_update
  }
  
  save(sample_mean_dist_MatrixRecord, file="sample_mean_dist_MatrixRecord.Rdata")
  colnames(sample_mean_dist_MatrixRecord)<-ks
  rownames(sample_mean_dist_MatrixRecord)<-sampleIDs
  write.table(sample_mean_dist_MatrixRecord, file="sample_mean_dist_MatrixRecord.xls" , sep="\t")
  

  #plot results and draw loess fitted line
  sample_mean_dist_MatrixRecord_avg<-colSums(sample_mean_dist_MatrixRecord)/nrow(sample_mean_dist_MatrixRecord)
  
  png('ElbowPlot.png')
  plot(ks, sample_mean_dist_MatrixRecord_avg, type="b", ylab="Average Within-Cluster Errors", xlab="Number of Clades", main="Elbow Point Analysis of Number of Clades")
  
  loess_result<-loess(sample_mean_dist_MatrixRecord_avg~ks)
  loess_ks <- seq(min(ks),max(ks), (max(ks) - min(ks))/1000)
  loess_SSE<-predict(loess_result, loess_ks)
  lines(loess_ks, loess_SSE, col='red', lwd=2)
  dev.off()
}
