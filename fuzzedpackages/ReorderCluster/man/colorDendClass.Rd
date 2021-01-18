\name{colorDendClass}
\alias{colorDendClass}

\title{
Makes the plot of the dendrogram, visualizing the class label information with different colors of dendrogram edges.
}
\description{
Makes the plot of the dendrogram, visualizing the class label information with different colors of dendrogram edges.
}
\usage{
colorDendClass(dend, class)
}

\arguments{
\item{dend}{the clustering dendrogram of the analyzed dataset.}
\item{class}{
a vector of length n, which stores the class labels of the dataset objects. 
}
}
\details{
Makes the plot of the dendrogram, visualizing the class label information with different colors of dendrogram edges.
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
dend=as.dendrogram(hc)

label=unique(class)

Rowcolor=rainbow(length(label))
rc=matrix(0,length(class),1)

for (j in 1:length(label))
{
		index=which(class==label[j])
		rc[index]=Rowcolor[j]
}

colorDendClass(dend,rc[hc$order])
}	

\keyword{ dendrogram }
\keyword{ clustering }
