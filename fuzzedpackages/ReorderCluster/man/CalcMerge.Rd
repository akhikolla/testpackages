\name{CalcMerge}
\alias{CalcMerge}

\title{
Forms the binary vector to mark the nodes with identical class labels.
}
\description{
This is the auxiliary function of the package. It generates the binary vector to mark the nodes with identical class labels.
}
\usage{
CalcMerge(hc,node,class)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hc}{
An object of class hclust which describes the tree produced by the clustering process.There are such components: merge, height, order, labels,call,method,dist.method.
}
  \item{node}{
a list of lists of length n-1 \code{\link{testBar}}, where each single list stores two vectors, consisting of elements of the left and right subtrees of the corresponding node in the merging matrix \code{hc}
}
 \item{class}{
a vector of length n, which stores the class labels of the dataset objects. 
}
}
\details{
This is the auxiliary function of the package. It generates the binary vector to mark the nodes with identical class labels. The node or subtree with the elements, having the same class labela are marked by value 1, otherwise 0. The output vector is used as input to function \code{\link{SubTree}}, which simplifies the hierarchical tree in order to increase the efficiency of the dynamic programming algorithm.
}
\value{
\item{flag}{
a binory vector of length n-1. 
}
}
\references{
Therese Biedl,Brona Brejova, Erik D,Demaine, Angele M.Hanmel and Tomas Vinar:Optimal Arrangement of Leaves in the Tree Representing Hierarchical Clustering of Gene Expression Data/Technical report 2001-14
}
\author{
Natalia Novoselova,Frank Klawonn,Junxi wang

Maintainer: Natalia Novoselova <novos65@mail.ru>
}

\seealso{
\code{\link{RearrangeJoseph}}, \code{\link{OrderingJoseph}}
}
\examples{
data(leukemia)
rownames(leukemia)=leukemia[,1]
leukemia=leukemia[,-1]
matr=leukemia[,-101]
class=leukemia[,101]

matr=as.matrix(matr)
dist=dist(matr)
hc <- hclust(dist)
node=testBar(hc)
flag=CalcMerge(hc,node,class)
}
\keyword{ clustering }
\keyword{ node }
