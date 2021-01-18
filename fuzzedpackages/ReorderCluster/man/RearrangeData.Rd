\encoding{UTF-8}
\name{RearrangeData}
\alias{RearrangeData}

\title{
Sample function to perform optimal reordering of the hierarchical tree according to class labels
}
\description{
Provides the main scheme to perform optimal reordering of the hierarchical tree according to class labels for any dataset
}
\usage{
RearrangeData(matr,class,clustering.method="complete",cDend=FALSE,dirpath="",cpp=TRUE)
}

\arguments{
\item{matr}{a dataset, a matrix of feature values for several cases.}

\item{class}{a character vector or numeric vector, containing the class labels of the dataset objects}
 \item{clustering.method}{
a method for performing the hierarchical clustering. The method names corresponds to the methods, that are used for the \code{hclust} function. The default value is "complete".
}
\item{cDend}{
A binary variable with values TRUE/FALSE enables or desable the visualiztion of the color dendrogram. Default is FALSE.
}
\item{dirpath}{
A path to the directory, where the heat map is saved. When the path is empry (default) the heat map is visualized on the display.
}
\item{cpp}{a binary variable, which allows to switch between dynamic programming algorithm, realized in R language \code{\link{OrderingJoseph}} and the same algorithm as the C++ function \code{\link{OrderingJosephC}}
}
}
\details{
The function visualizes the heat map for the results of the hiererchical clustering of the input dataset, than the call to the function \code{\link{RearrangeJoseph}} enables the reordering of the dendrogram according to the supplied class labels. At the end the resultant dendrogram and the heat map are visualized.
}

\value{
\item{order}{a numeric vector of the reordered leaves of the debdrogram.
}
\item{A}{a numerical nxn matrix, which stores the values of the evaluation function for the subtrees. Each cell (i,j) of the matrix A stores the best ordering (according to the evaluation function value) of the subtree starting with element i and ending with element j.
}
\item{class}{a vector of labels of the dataset objects, reordered according to the vector \code{order}.
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
res=RearrangeData(matr,class)
}  

\keyword{ clustering }
\keyword{ reordering }
