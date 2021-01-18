
#' Creates network alignments using network constructed from subpopulations after PAC
#' 
#' @param sampleIDs sampleID vector
#' @param num_PACSupop number of subpopulations learned in PAC step for each sample
#' @param smallSubpopCutoff Population size cutoff for subpopulations in clade calculation. The small subpopulations will be considered in the refinement step.
#' @param k_clades number of clades to output before refinement
#' @return clades_network_only the clades constructed without small subpopulations (by cutoff) using mutual information network alignments
#' @export

MAN<-function(sampleIDs, num_PACSupop, smallSubpopCutoff, k_clades){
  
  Networks<-NULL
  data_agg<-NULL
  
  mainDir<-getwd()
  #create network dissimilarity matrix
  network_list<-list(0)
  n=1
  for(i in 1:length(sampleIDs)){
    subDir <- paste0("/", sampleIDs[i], "_SampleSubpopulationMatrixNetworks")
    setwd(paste0(mainDir, subDir))
    
    for(j in 1:num_PACSupop){
      load(paste0(sampleIDs[i], "_", j, "_subpop_Network.Rdata"))
      network_list[[n]]<-Networks
      names(network_list)[[n]]<-paste(sampleIDs[i], "_", j, sep="") 
      n=n+1
    }
    
    setwd(mainDir)
    
  }
  
  #dissimilarity calculation using Jaccard distance
  all_Networks_dissim<-matrix(0,length(network_list),length(network_list))
  colnames(all_Networks_dissim)<-names(network_list)
  rownames(all_Networks_dissim)<-names(network_list)
  
  for(q in 1:length(network_list)){
    for(w in 1:length(network_list))
      all_Networks_dissim[q,w]<-JaccardSM(network_list[[q]], network_list[[w]])
  }
  
  #alignments of networks without very small subpopulations
  all_data_agg<-matrix(0,0,0)
  for(r in 1:length(sampleIDs)){
    load(paste0(sampleIDs[r],"_data_agg.Rdata"))
    all_data_agg<-rbind(all_data_agg, data_agg)
  }
  
  subpopNames<-paste0(all_data_agg[,2], "_", all_data_agg[,1])
  all_data_agg<-all_data_agg[,-c(1,2)]
  rownames(all_data_agg)<-subpopNames
  
  subpopSizeFilter<-all_data_agg[,ncol(all_data_agg)]>=smallSubpopCutoff
  xN<-all_Networks_dissim[subpopSizeFilter,subpopSizeFilter]
  
  xN<-1-xN
  xN<-as.dist(xN)
  hcN<-hclust(xN, "ave")
  
  png("hcN_clusterDendrogram_subset.png")
  plot(hcN)
  dev.off()
  
  groupA<-cutree(hcN, k = k_clades)
  clades_network_only<-split(groupA, groupA)
  return(clades_network_only)
}