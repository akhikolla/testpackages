#' @importFrom dplyr %>% group_by summarise_all funs
NULL

#' Creates the matrix that can be easily plotted with a heatmap function available in an R package
#' 
#' @param aggregateMatrix_withAnnotation the annotated clade matrix
#' @return the heatmap input matrix
#' @export


heatmapInput<-function(aggregateMatrix_withAnnotation){
  
  ClusterID<-NULL
  SampleID<-NULL
  
  data_agg<-aggregateMatrix_withAnnotation
  cluster_index <- levels(data_agg[,2])
  sample_names <- unique(data_agg[, 3])
  cluster_index_pad <- rep(cluster_index, length(sample_names))
  sample_names_pad <- rep(sample_names, each = length(cluster_index))
  
  Zero_Pad <- matrix(0, length(sample_names_pad), 1)
  ZeroPadding <- data.frame(cluster_index_pad, sample_names_pad, Zero_Pad, stringsAsFactors = FALSE)
  dat_count_samples <- data_agg[, c(2, 3, ncol(data_agg))]
  
  
  colnames(ZeroPadding) <- colnames(dat_count_samples)
  padded_dat_count_samples <- rbind.data.frame(ZeroPadding,dat_count_samples, stringsAsFactors = FALSE)
  data_agg_intermediate <- padded_dat_count_samples %>% group_by(ClusterID,  SampleID) %>% summarise_all(funs(sum))
  data_agg_intermediate <- data.frame(data_agg_intermediate, stringsAsFactors = FALSE)
  
  cladeProportionMatrix<-matrix(0,length(sample_names),length(cluster_index))
  cladeProportionMatrix<-as.data.frame(cladeProportionMatrix)
  rownames(cladeProportionMatrix)<-sample_names
  for(i in 1:length(cluster_index)){
    data_agg2 <- data_agg_intermediate[order(data_agg_intermediate$SampleID),]
    clade_count_summary<-data_agg2[which(data_agg2[,1]==cluster_index[i]),]
    clade_proportion_summary<-clade_count_summary[,3]/sum(clade_count_summary[,3])
    cladeProportionMatrix[,i]<-clade_proportion_summary
    colnames(cladeProportionMatrix)[i]<-cluster_index[i]
  }
  
  return(cladeProportionMatrix)
  
}
