\encoding{UTF-8}
\name{ReorderCluster-package}
\alias{ReorderCluster-package}
\alias{ReorderCluster}
\docType{package}
\title{
optimal reordering of the hierarchical tree according to class labels
}
\description{
The package includes the functions to make the near optimal reordering of the leaves of dendrogram  according to the class labels, preserving its initial structure. The reordering algorithm tries to group instances from the same class together. With such an optimised dendrogram it is easier to interpret the interrelation between the available class labels and clusters. Moreover, it is possible to see whether one can predict the class for unlabeled instances based on the distance that was used for clustering.
}
\details{
\tabular{ll}{
Package: \tab ReorderCluster\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-07-21\cr
License: \tab GPL (>= 3)\cr
}
The basic unit of ReorderCluster package is the function \code{\link{RearrangeJoseph}}, which makes the initialization of all the auxiliary matrices and the sequence of the necessary functions' call to perform the leaf reordering according to class labels. The function \code{\link{RearrangeJoseph}} implements the preprocessing of the merging matrix of the hierarchical clustering by calling the auxiliary function \code{\link{testBar}} to form the node (subtree) structure using the merging matrix, the auxiliary functions \code{\link{CalcMerge}} to mark the nodes with identical class labels, the auxiliary functions \code{\link{SubTree}} to simplify the initial hierarchical tree by constructing the modified tree, where the subtrees with equal class labels are merged into the single element. The function \code{\link{RearrangeJoseph}} calls the function \code{\link{OrderingJoseph}} to execute the dynamic programming algorithm, which calculates the values of the evaluation function for each subtree of the dendrogram and thereafter finds the optimized sequence of leaves using the function \code{\link{funMerge}}. The evaluation function is calculated on the basis of the existing class labels.
}

\section{Functions}{
  \tabular{ll}{
    \code{\link{RearrangeJoseph}} \tab Makes the calls to sequence of functions to perform the leaf reordering according to class labels.\cr
    \code{\link{OrderingJoseph}} \tab Calculates the evaluation function for each subtree using the dynamic programming approach.\cr
    \code{\link{OrderingJosephC}} \tab Calculates the evaluation function for each subtree using the dynamic programming approach calling the function in C++ language.\cr
    \code{\link{RearrangeData}} \tab Provides the main scheme to perform optimal reordering of the hierarchical tree according to class labels for any dataset.\cr
    \code{\link{funMerge}} \tab Recovers the optimal sequence of leaves of the hierarchical tree using the backtracking scheme.\cr
    \code{\link{CalcMerge}} \tab Forms the auxiliary vectors to mark the nodes with identical class labels.\cr
    \code{\link{SubTree}} \tab Simplifies the initial hierarchical tree by reducing the number of nodes.\cr
    \code{\link{testBar}} \tab For each inner node of the dendrogram forms two vectors with left and write subtrees' elements.\cr
    \code{\link{colorDendClass}} \tab Makes the plot of the dendrogram, visualizing the class label information with different colors of dendrogram edges.\cr
    \code{\link{testData1}} \tab Generates the simulated dataset with 400 genes and 100 experiments with 3 class labels.\cr
    \code{\link{testData2}} \tab Generates the simulated dataset with 90 genes and 90 experiments with 3 class labels.\cr
  }
}
\section{Installing and using}{
  To install this package, make sure you are connected to the internet and issue the following command in the R prompt:
  
  \preformatted{
    install.packages("ReorderCluster")
  }
  
  To load the package in R:
 
  \preformatted{
    library(ReorderCluster)
  }
}
\author{
Natalia Novoselova,Frank Klawonn,Junxi wang

Maintainer: Natalia Novoselova <novos65@mail.ru>
}
\references{
Z. Bar-Joseph, D.K. Gifford, and T.S. Jaakkola. Fast optimal leaf ordering for hierarchical clustering. Bioinformatics, 17:22-29, 2001.
Therese Biedl,Brona Brejova, Erik D,Demaine, Angele M.Hanmel and Tomas Vinar:Optimal Arrangement of Leaves in the Tree Representing Hierarchical Clustering of Gene Expression Data/Technical report 2001-14.
}
\keyword{ package }
\keyword{clustering}
\keyword{reorder}
\keyword{optimal}

\seealso{
CRAN packages \pkg{cba}, \pkg{seriation}
includes the function for optimal leaf reordering of the dendrogram with respect to the minimizing the sum of the distances along the (Hamiltonian) path connecting the leaves in the given order.
}
\examples{
  data(leukemia)
  cpp=FALSE
  rownames(leukemia)=leukemia[,1]
  leukemia=leukemia[,-1]
  matr=leukemia[,-101]
  class=leukemia[,101]
  
  matr=as.matrix(matr)
  dd=dim(matr)

	label=unique(class)
	
	Rowcolor=rainbow(length(label))
	rc=matrix(0,length(class),1)
	ds=matrix(0,length(class),1)

	for (j in 1:length(label))
  	{
		index=which(class==label[j])
		rc[index]=Rowcolor[j]
	}
	
	dist=dist(matr)
	hc <- hclust(dist)
	
	dend=as.dendrogram(hc)
  
	my_palette <- colorRampPalette(c("green", "black", "red"))(n = 399)
	hv <- heatmap.2(matr,Rowv=dend,scale = "none",Colv=NA,
	                col=my_palette, RowSideColors = rc,trace="none",dendrogram="row")
  legend("topright",legend=label,col=Rowcolor,pch=15,cex=0.8)

  res=RearrangeJoseph(hc,as.matrix(dist),class,cpp) 
  
	hcl=res$hcl

  dend=as.dendrogram(hcl)
  
	hv <- heatmap.2(matr,Rowv=dend,scale = "none",Colv=NA,
	                col=my_palette, RowSideColors = rc,trace="none",dendrogram="row")
  legend("topright",legend=label,col=Rowcolor,pch=15,cex=0.8)
}
