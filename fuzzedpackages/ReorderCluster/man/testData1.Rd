\name{testData1}
\alias{testData1}

\title{
Simulates the dataset for analysis.
}
\description{
Generates the simulated dataset with 400 genes and 100 experiments with 3 class labels.
}
\usage{
testData1()
}


\details{
Generates the simulated dataset with 400 genes and 100 experiments with 3 class labels. The last column stores the class labels for dataset objects.
}
\value{
\item{matr}{a dataset, a matrix of feature values for several cases. The additional last column stores the values of class labels.}
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
matr=testData1()
}
\keyword{ dataset }
\keyword{ clustering }
