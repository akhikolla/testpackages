\name{dataSimHuman}
\docType{data}
\alias{dataSimHuman}
\title{Simulated cell signalling pathway data}
\description{
 Data generated from the known regulatory network of human cell signalling data. 
}
\usage{data(SimHumanPw)}
\format{
  \code{dataSimHuman} contains the following objects:
  \describe{
    \item{\code{Obs}}{
      List of 100 matrices (100x11) each of which with 100 observations simulated from the known regulatory network.
    }
    \item{\code{Perms}}{
      List of 5 matrices (1x11) each of which with a permutation of the nodes.
    }
    \item{\code{TDag}}{
      Matrix (11x11) with the adjacency matrix of the known regulatory network.
    }
  }
}
\source{ 
D. Altomare, G. Consonni and L. La Rocca (2012). Objective Bayesian search of Gaussian directed acyclic graphical models for ordered variables with non-local priors. \emph{Article submitted to Biometric Methodology}.
}
\references{
Sachs, K., Perez, O., Pe'er, D., Lauffenburger, D., and Nolan, G. (2003). Casual protein-
signaling networks derived from multiparameter single-cell data. \emph{Science} 308, 504-6.

Shojaie, A. and Michailidis, G. (2010). Penalized likelihood methods for estimation of sparse
high-dimensional directed acyclic graphs. \emph{Biometrika} 97, 519-538.
}
\keyword{datasets}
