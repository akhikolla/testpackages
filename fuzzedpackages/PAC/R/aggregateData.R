#' @importFrom dplyr %>% group_by summarise_all funs
NULL

#' Aggregates results from the clustering and merging step.
#' 
#' @param dataInput Data matrix, with first column being SampleID.
#' @param labelsInput cluster labels from PAC.
#' @return The aggregated data of \code{dataInput}, with average signal levels for all clusters and sample combinations.
#' @export
#' @examples 
#' n = 5e3                       # number of observations
#' p = 1                         # number of dimensions
#' K = 3                         # number of clusters
#' w = rep(1,K)/K                # component weights
#' mu <- c(0,2,4)                # component means
#' sd <- rep(1,K)/K              # component standard deviations
#' g <- sample(1:K,prob=w,size=n,replace=TRUE)   # ground truth for clustering
#' X <- as.matrix(rnorm(n=n,mean=mu[g],sd=sd[g]))
#' y <- PAC(X, K)
#' X2<-as.matrix(rnorm(n=n,mean=mu[g],sd=sd[g]))
#' y2<-PAC(X2,K)
#' X<-cbind("Sample1", as.data.frame(X)); colnames(X)<-c("SampleID", "Value")
#' X2<-cbind("Sample2", as.data.frame(X2)); colnames(X2)<-c("SampleID", "Value")
#' aggregateData(rbind(X,X2),c(y,y2))


aggregateData<-function(dataInput, labelsInput){
  ClusterID<-NULL
  SampleID<-NULL
  
  data<-cbind(labelsInput, dataInput)
  colnames(data)[1]<-"ClusterID"
  colnames(data)[2]<-"SampleID"
  data<-data.frame(data, stringsAsFactors = FALSE)
  data$count<-1
  data_agg_intermediate<-data %>% group_by(ClusterID, SampleID) %>% summarise_all(funs(sum))
  data_agg_intermediate<-data.frame(data_agg_intermediate)
  data_agg<- data_agg_intermediate[order(data_agg_intermediate$SampleID),]
  data_agg[,-c(1,2,ncol(data_agg))]<-data_agg[,-c(1,2,ncol(data_agg))]/data_agg$count
  
  return(data_agg)
  
}