\encoding{latin1}
\name{centerX}
\alias{centerX}

\title{Cluster-center a design matrix around cluster means}

\description{Given a design matrix and a cluster identification
  variable. A new design matrix is returned, where each column
  is centered around its cluster mean. Using C++ code for speed.}
\usage{
centerX(x, id)
}
\arguments{
  \item{x}{A design matrix}
  \item{id}{A cluster identifying variable.}
}
\details{
  \code{centerX} is written as a help function with the aim of being fast.}
\value{ A matrix with cluster centered columns.}
\author{
Johan Zetterqvist, Arvid \enc{Sjölander}{Sjolander}
}
\keyword{cluster}
