\name{MU_Fast}
\alias{MU_Fast}
\title{
Fast Implementation of First Order Moment of Generalized Hamming Distance Statistic
}
\description{
Provides a very fast implementation for obtaining the asymptotic mean for GHD statistic
}
\usage{
MU_Fast(A, B)
}
\arguments{
  \item{A}{
   Topological matrix obtained from graph g_A.
}
  \item{B}{
   Topological matrix obtained from graph g_B.
}
}
\value{
  Returns the asymptotic mean value for GHD statistic between permutations of topological matrices of graphs g_A and graph g_B.
}
\author{
Raghvendra Mall <rmall@hbku.edu.qa>
}

\seealso{
\code{\link{GHD_Fast}}, \code{\link{STD_Fast}}
}
