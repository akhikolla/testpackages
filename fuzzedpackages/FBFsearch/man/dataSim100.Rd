\name{dataSim100}
\docType{data}
\alias{dataSim100}
\title{DAG model with 100 nodes and 100 edges}
\description{
\code{dataSim100} is a \code{list} with the adjacency matrix of a randomly generated DAG with 100 nodes and 100 edges, 10 samples generated from the DAG and 5 permutations of the nodes.}
\usage{data(SimDag100)}
\format{
  \code{dataSim100} contains the following objects:
  \describe{
    \item{\code{Obs}}{
      List of 10 matrices (100x100) each of which with 100 observations generated from the DAG.
    }
    \item{\code{Perms}}{
      List of 5 matrices (1x100) each of which with a permutation of the nodes.
    }
    \item{\code{TDag}}{
      Matrix (100x100) with the adjacency matrix of the DAG.
    }
  }
}
\source{ D. Altomare, G. Consonni and L. La Rocca (2012). Objective Bayesian search of Gaussian directed acyclic graphical models for ordered variables with non-local priors. \emph{Article submitted to Biometric Methodology}.
}
\references{Shojaie, A. and Michailidis, G. (2010). Penalized likelihood methods for estimation of sparse
high-dimensional directed acyclic graphs. \emph{Biometrika} 97, 519-538.
}
\keyword{datasets}
