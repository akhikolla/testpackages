\name{groupcount}
\alias{groupcount}
\title{
Generate 2D count table for two integer-valued vectors.
}
\description{
Similar to 'table' function, this function builts a contingency table of the counts at each combination of all possible values from two integer-valued input vectors.
}
\usage{
groupcount(g1, g2, n1, n2)
}
\arguments{
  \item{g1}{
  The first integer-valued input vector. The max value in g1 is n1.
}
  \item{g2}{
  The second integer-valued input vector. The max value in g1 is n2.
}
  \item{n1}{
  The maximum value in g1.
}
  \item{n2}{
  The maximum value in g2.
}
}
\details{
This is implemented as an utility function to build a 2D histogram count table. For efficiency, it does not check if the maximum values in input vectors exceed the maximum values specified.
}
\value{
The count table.
}

\author{
Quanli Wang
}

\examples{
n1 <- 20
n2 <- 10
g1 <- sample.int(n1,1000, replace = TRUE)
g2 <- sample.int(n2,1000, replace = TRUE)
counts <- groupcount(g1,g2,n1,n2)
}
