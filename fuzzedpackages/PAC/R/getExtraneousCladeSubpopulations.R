
#' Calculates subpopulations in clades (with two or more subpopulations) that are too far away from other subpopulations (within the same clade) on the constellation plot; these far away subpopulations should be pruned away from the original clades.
#' @param tsneResults t-SNE output of clade centroids' embedding.
#' @param pacman_results PAC-MAN analysis result matrix that contains network annotation, clade IDs and mean (centroid) clade expression levels.
#' @param threshold_multiplier how many times the threshold ( (a) spread from center of clade for clades with three or more sample subpopulations and (b) distance from each subpopulation centroid for clades with exactly two subpopulations).     
#' @param max_threshold the maximum distance (on t-SNE plane) allowed for sample subpopulations to be categorized into the same clade.
#' 
#' @return Returns clade subpopulations to be pruned.
#' 
#' @export



getExtraneousCladeSubpopulations<-function(tsneResults, pacman_results,threshold_multiplier,max_threshold){
  
  ClusterID<-NULL
  
  clade_subset_to_change<-matrix(NA, nrow=0, ncol=2)
  colnames(clade_subset_to_change)<-c("ClusterID", "SampleID")
  
  combined_tsne_pacman_results<-cbind(as.data.frame(tsneResults), pacman_results[,c(2,3)])
  
  #clades that have exactly 2 subpopulations
  recurrentClades_2<-names(which(table(combined_tsne_pacman_results[,3]) == 2))
  combined_tsne_pacman_results_2<-combined_tsne_pacman_results[(combined_tsne_pacman_results[,3] %in% recurrentClades_2),]
  
  clade_subset<-vector()
  
  for(u in 1:length(recurrentClades_2)){
    clade_subset<-subset(combined_tsne_pacman_results_2, ClusterID == recurrentClades_2[u])
    
    clade_subset_info<-clade_subset[,c(3,4)]
    clade_subset<-clade_subset[,c(1,2)]
    
    #distance between two points
    clade_dist<-sqrt(sum(clade_subset[1,]-clade_subset[2,])^2)
    threshold<-threshold_multiplier*getAverageSpreadOf2SubpopClades(tsneResults, pacman_results)
    
    threshold<-min(threshold, max_threshold)
    
    #cladeRange<-c(clade_avgDist-threshold, clade_avgDist+threshold)
    #outlier<-euclidean_dist_fromCentroid < cladeRange[1]
    
    if(clade_dist > threshold){
      clade_subset_to_change<-rbind(clade_subset_to_change,clade_subset_info[1,]) #prune away 1 of the 2 subpopulations
    }
    
  }
  
  
  #clades that have more than two subpopulations
  recurrentClades_3<-names(which(table(combined_tsne_pacman_results[,3]) >= 3))
  combined_tsne_pacman_results_3<-combined_tsne_pacman_results[(combined_tsne_pacman_results[,3] %in% recurrentClades_3),]
  
  clade_subset<-vector()
  for(i in 1:length(recurrentClades_3)){
    clade_subset<-subset(combined_tsne_pacman_results_3, ClusterID == recurrentClades_3[i])
    
    clade_subset_info<-clade_subset[,c(3,4)]
    clade_subset<-clade_subset[,c(1,2)]
    
    matrix_centroid<-colSums(clade_subset)/nrow(clade_subset)
    
    euclidean_dist_fromCentroid<-vector()
    for(j in 1:nrow(clade_subset)){
      euclidean_dist_fromCentroid[j]<-sqrt(sum(clade_subset[j,]-matrix_centroid)^2)
    }
    
    clade_avgDist<-mean(euclidean_dist_fromCentroid)
    threshold<-threshold_multiplier*(clade_avgDist)
    
    threshold<-min(threshold, max_threshold)
    
    cladeRange<-c(clade_avgDist-threshold, clade_avgDist+threshold)
    outliersLow<-euclidean_dist_fromCentroid < cladeRange[1]
    outliersHigh<-euclidean_dist_fromCentroid > cladeRange[2]
    
    if(sum(outliersLow) + sum(outliersHigh) > 0){
      clade_subset_to_change<-rbind(clade_subset_to_change,clade_subset_info[which(as.numeric(outliersHigh)==1),])
      clade_subset_to_change<-rbind(clade_subset_to_change,clade_subset_info[which(as.numeric(outliersLow)==1),])
    }
    
  }
  
  return(clade_subset_to_change)
}


