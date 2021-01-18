% File src/library/base/man/load.Rd
\name{CollapseLabels}
\alias{CollapseLabels}
\title{CollapseLabels}
\description{
  Performs a permutation on the group labels so that, if \code{K} non-empty groups are present, the labels used are exactly \{\code{1, \dots, K}\}.
}
\usage{
CollapseLabels(decision)
}
\arguments{
  \item{decision}{A vector identifying a partition. The entries must be positive integers.}
}
\details{
The labels of the clustering provided \code{{1,\dots,G}} are mapped into \code{{1,\dots,K}} (\code{K is less or equal than G}) based on their order of appearance in the vector provided.
}
\examples{
set.seed(123)
decision <- sample(1:50, 50, TRUE)
CollapseLabels(decision = decision)
}



