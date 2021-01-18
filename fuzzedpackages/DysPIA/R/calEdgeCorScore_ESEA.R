#' @title calEdgeCorScore_ESE
#' @description Calculates differential Mutual information.
#'
#' @param dataset Matrix of gene expression values (rownames are genes, columnnames are samples).
#' @param class.labels Vector of binary labels.
#' @param controlcharacter Charactor of control in the class labels.
#' @param casecharacter Charactor of case in the class labels.
#' @param background Matrix of the edges' background.
#' @return A vector of the aberrant correlation in phenotype P based on mutual information (MI) for each edge.
#' @importFrom parmigene knnmi
#' @importFrom stats na.omit var
#' @import DysPIAData
#' @examples
#' data(gene_expression_p53, class.labels_p53,sample_background)
#' ESEAscore_p53<-calEdgeCorScore_ESEA(gene_expression_p53, class.labels_p53,
#'  "WT", "MUT", sample_background)
#' 
calEdgeCorScore_ESEA <- function(dataset, class.labels, controlcharacter, casecharacter, background) {

  controlloca<-which(class.labels==controlcharacter)
  caseloca<-which(class.labels==casecharacter)
  contrSet<-dataset[,controlloca]
  datasetloca<-c(controlloca,caseloca)
  dataset<-dataset[,datasetloca]

  data_var<-apply(dataset,1,var)*apply(contrSet,1,var)
  if (sum(data_var==0)>0){dataset<-dataset[-which(data_var==0),]}
  
location<-matrix(0,length(background[,1]),2)
location[,1]<-match(background[,1],rownames(dataset))
location[,2]<-match(background[,2],rownames(dataset))
location<-na.omit(location)

dataset.1<-dataset[location[,1],]
dataset.2<-dataset[location[,2],]
Cexpress.1<-dataset[location[,1],1:ncol(contrSet)]
Cexpress.2<-dataset[location[,2],1:ncol(contrSet)]
EdgeCorScore<-c()

for(i in 1:length(location[,1])){
  EdgeCorScore[i]<-knnmi(dataset.1[i,],dataset.2[i,])-knnmi(Cexpress.1[i,],Cexpress.2[i,])
}
EdgeID<-matrix(0,length(location[,1]),2)
EdgeID[,1]<-row.names(dataset.1)
EdgeID[,2]<-row.names(dataset.2)
colnames(EdgeID)<-c("Edge1","Edge2")
Paste<-function(x,c1,c2) paste(x[c1],x[c2],sep="|")
names(EdgeCorScore)<-apply(EdgeID,1,Paste,c1="Edge1",c2="Edge2")
return(EdgeCorScore)
}
