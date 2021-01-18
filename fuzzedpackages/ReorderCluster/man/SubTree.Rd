\name{SubTree}
\alias{SubTree}

\title{
Simplifies the initial hierarchical tree by reducing the number of nodes. Constructs the new merging matrix with some inner nodes substituted by one element from the coressponding subtree.
}
\description{
This is the auxiliary function of the package. It creates the four-column matrix for reordering, that contains the results of merging the subtrees of the dendrogram with identical class labels into single elements.The result is the simplifies tree, which allows to spare the computational time and to increase the efficiency of the reordering algorithm \code{\link{OrderingJoseph}}.
}
\usage{
SubTree(ind,flag,node,hc,A,r,coef)
}

\arguments{
  \item{ind}{
the value of the last row of the merge matrix, which is equal to n-1, where n is the number of data objects.
}
\item{flag}{
a binory vector of length n-1. 
}
  \item{node}{
a list of lists of length n-1 \code{\link{testBar}}, where each single list stores two vectors, consisting of elements of the left and right subtrees of the corresponding node in the merging matrix \code{hc}
}
\item{hc}{
An object of class hclust which describes the tree produced by the clustering process.There are such components: merge, height, order, labels,call,method,dist.method.
}
\item{A}{
a numerical nxn matrix, which stores the values of the evaluation function fir the subtrees. Each cell (i,j) of the matrix A stores the best ordering (according to the evaluation function value) of the subtree starting with element i and ending with element j.
}
\item{r}{
a numerical nxn matrix, each cell stores the number of the elements with the same class label, starting from the leftmost or the rightmost element of the subtree.
}
\item{coef}{
a parameter used in the expression for the calculation of the evaluation function. It reinforces the influence of longer sequences. Its values are in the interval ]1,2].
}
}
\details{
This function modifies the merging matrix of \code{hc}, by merging the subtrees, which have the same elements' class labels. As a result the n-1 by 2 matrix is formed with the columns coressponding to nodes of the simplified hierarchical tree. Using this merging matrix allows to rise the time efficiency of the reordering algorithm (see \code{\link{OrderingJoseph}} and \code{\link{RearrangeJoseph}}).
}

\value{
\item{hc}{
An object of class hclust which include the modified merging matrix hc$merge, with the subtrees with identical class label merged into a single element (the leftmost element of the corresponding subtree)
}
\item{node}{
the modified input structure \code{node} (a list of lists of length n-1), which corresponds to the modified \code{hc}
}
\item{A}{
the modified nxn matrix \code{A}, which corresponds to the modified \code{hc}.
}
\item{r}{
a modified nxn matrix \code{r}, which corresponds to the modified \code{hc}.
}
}
\references{
Z. Bar-Joseph, D.K. Gifford, and T.S. Jaakkola. Fast optimal leaf ordering for hierarchical clustering. Bioinformatics, 17:22-29, 2001.
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

coef=1.5
num=dim(hc$merge)[1]
A=array(1,c(num+1,num+1))
r=array(1,c(num+1,num+1))

ind=num
node=testBar(hc)
flag=CalcMerge(hc,node,class)
out=SubTree(ind,flag,node,hc,A,r,coef)
}
\keyword{ clustering }
\keyword{ node }
