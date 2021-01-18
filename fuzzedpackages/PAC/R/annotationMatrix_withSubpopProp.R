#' Adds subpopulation proportion for the annotation matrix for the clades
#' 
#' @param aggregateMatrix_withAnnotation the annotated clade matrix
#' @return Annotated clade matrix with subpopulation proportions
#' @export


annotationMatrix_withSubpopProp<-function(aggregateMatrix_withAnnotation){
  annotationMatrix_prop<-matrix(0,0,0)
  for(i in 1:length(unique(aggregateMatrix_withAnnotation[,3]))){
    sampleID_i<-unique(aggregateMatrix_withAnnotation[,3])[i]
    sampleMeans<-aggregateMatrix_withAnnotation[which(aggregateMatrix_withAnnotation[,3]==sampleID_i),]
    counts<-sampleMeans[,ncol(sampleMeans)]
    subpop_proportion<-counts/sum(counts)
    subpop_proportion<-signif(subpop_proportion*100, digits=4)
    annotationMatrix_prop<-rbind(annotationMatrix_prop,cbind(sampleMeans,subpop_proportion))
  }
  
  write.table(annotationMatrix_prop, file="annotationMatrix_withSubpopProp.xls", sep="\t", row.names=FALSE)
  return(annotationMatrix_prop)
}
