\name{testBar}
\alias{testBar}

\title{
For each node (subtree) of the hierarchical tree forms two vectors, consisting of elements of its left and write subtrees.
}
\description{
This is the auxiliary function of the package. It makes the lists of two vector, consisting of elements of left and write subtrees for each node of the hierarchical tree. These lists are used in reordering algorithm \code{\link{OrderingJoseph}}, so that it makes easy to define the first(leftmost) and last(rightmost) element of each subtree.
}
\usage{
testBar(hc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hc}{
An object of class hclust which describes the tree produced by the clustering process.There are such components: merge, height, order, labels,call,method,dist.method.
}
}
\details{
For each node or subtree of the hierarchical clustering result the function forms the list of two vector, consisting of elements of its left and write subtrees. These lists are used in reordering algorithm \code{\link{OrderingJoseph}}, so that it makes easy to define the first(leftmost) and last(rightmost) element of each subtree, and also to determine the intermediate elements in the process of searching the best ordering for each particular subtree.
}
\value{
\item{node}{a list of lists of length n-1, where each single list stores two vectors, consisting of elements of the left and right subtrees of the corresponding node in the merging matrix \code{hc}}
}
\references{
Z. Bar-Joseph, D.K. Gifford, and T.S. Jaakkola. Fast optimal leaf ordering for hierarchical clustering. Bioinformatics, 17:22-29, 2001.
}
\author{
Natalia Novoselova,Frank Klawonn,Junxi wang

Maintainer: Natalia Novoselova <novos65@mail.ru>
}

\seealso{
\code{\link{CalcMerge}}, \code{\link{SubTree}}
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
}
\keyword{ clustering }
\keyword{ node }
