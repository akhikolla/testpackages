\name{CollapseLabels}
\alias{CollapseLabels}
\title{CollapseLabels}
\description{
  Given a matrix of allocations, it performs a permutation on the group labels so that, if \code{K} non-empty groups are present, the labels used are exactly \{\code{1, \dots, K}\}. Note that the value \code{z[t,i]=0} is reserved to the inactive nodes: this transformation leaves such allocations unchanged.
}
\usage{
CollapseLabels(allocations)
}
\arguments{
  \item{allocations}{A matrix whose rows identify partitions of the same elements. The entries must be positive integers. Entries equal to zero are interpreted as missing values and hence are left unchanged.}
}
\details{
The labels \code{{1,\dots,G}} of the clustering provided are mapped into \code{{1,\dots,K}} (\code{K is less or equal than G}) based on their order of appearance in the matrix provided.
}
\examples{
set.seed(12345)
allocations <- matrix(sample(0:25, 25, TRUE), 5, 5)
CollapseLabels(allocations = allocations)
}


