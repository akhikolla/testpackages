\name{dataSim6}
\docType{data}
\alias{dataSim6}
\title{DAG model with 6 nodes and 5 edges}
\description{
\code{dataSim6} is a \code{list} with the adjacency matrix of a randomly generated DAG with 6 nodes and 5 edges and 100 correlation matrices generated from the DAG.
}
\usage{data(SimDag6)}
\format{
  \code{dataSim6} contains the following objects:
  \describe{
    \item{\code{Corr}}{
      List of 100 matrices (6x6) each of which with a correlation matrix generated from the DAG.
    }
    \item{\code{TDag}}{
      Matrix (6x6) with the adjacency matrix of the DAG.
    }
  }
}
\references{
D. Altomare, G. Consonni and L. La Rocca (2012). Objective Bayesian search of Gaussian directed acyclic graphical models for ordered variables with non-local priors. \emph{Article submitted to Biometric Methodology}.
}
\keyword{datasets}
