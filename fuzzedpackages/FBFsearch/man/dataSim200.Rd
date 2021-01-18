\name{dataSim200}
\docType{data}
\alias{dataSim200}
\title{DAG model with 200 nodes and 100 edges}
\description{
 \code{dataSim200} is a \code{list} with the adjacency matrix of a randomly generated DAG with 200 nodes and 100 edges, 10 samples generated from the DAG and 5 permutations of the nodes.}
\usage{
data(SimDag200)
}
\format{
  \code{dataSim200} contains the following objects:
  \describe{
    \item{\code{Obs}}{
      List of 10 matrices (100x200) each of which with 100 observations simulated from the DAG.
    }
    \item{\code{Perms}}{
      List of 5 matrices (1x200) each of which with a permutation of the nodes.
    }
    \item{\code{TDag}}{
      Matrix (200x200) with the adjacency matrix of the DAG.
    }
  }
}
\source{ D. Altomare, G. Consonni and L. La Rocca (2012). Objective Bayesian search of gaussian directed acyclic graphical models for ordered variables with non-local priors. \emph{Article submitted to Biometric Methodology}.
}
\references{Shojaie, A. and Michailidis, G. (2010). Penalized likelihood methods for estimation of sparse
high-dimensional directed acyclic graphs. \emph{Biometrika} 97, 519-538.
}
\keyword{datasets}
