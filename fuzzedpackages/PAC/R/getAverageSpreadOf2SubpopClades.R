
#' Calculate the (global) average spread of subpopulations in clades with 2 subpopulations on the constellation plot.
#' @param tsneResults t-SNE output of clade centroids' embedding.
#' @param pacman_results PAC-MAN analysis result matrix that contains network annotation, clade IDs and mean (centroid) clade expression levels.
#' 
#' @return Returns global average of 2-subpopulation clade spread on the constellation plot.



getAverageSpreadOf2SubpopClades<-function(tsneResults, pacman_results){
  ClusterID<-NULL
  
  l<-cbind(as.data.frame(tsneResults), pacman_results[,c(2,3)])
  recurrentClades<-names(which(table(l[,3]) == 2))
  l2<-l[(l[,3] %in% recurrentClades),]
  
  clade_dist_record<-vector()
  
  for(i in 1:length(recurrentClades)){
    
    clade_subset<-subset(l2, ClusterID == recurrentClades[i])
    
    clade_subset_info<-clade_subset[,c(3,4)]
    clade_subset<-clade_subset[,c(1,2)]
    
    #distance between two points
    clade_dist<-sqrt(sum(clade_subset[1,]-clade_subset[2,])^2)
    clade_dist_record[i]<-clade_dist
  }
  
  avgDist<-mean(clade_dist_record)
  return(avgDist)
}
