\name{dataHuman}
\docType{data}
\alias{dataHuman}
\title{Cell signalling pathway data}
\description{
 Data on a set of flow cytometry experiments on signaling networks of human immune system cells. The dataset includes p=11 proteins and n=7466 samples.
}
\usage{data(HumanPw)}
\format{
  \code{dataHuman} contains the following objects:
  \describe{
    \item{\code{Obs}}{
      Matrix (7466x11) with the observations.
    }
    \item{\code{Perms}}{
      List of 5 matrices (1x11) each of which with a permutation of the nodes.
    }
    \item{\code{TDag}}{
      Matrix (11x11) with the adjacency matrix of the known regulatory network.
    }
  }
}
\source{ Sachs, K., Perez, O., Pe'er, D., Lauffenburger, D., and Nolan, G. (2003). Casual protein-
signaling networks derived from multiparameter single-cell data. \emph{Science} 308, 504-6.
}

\references{
D. Altomare, G. Consonni and L. La Rocca (2012). Objective Bayesian search of Gaussian directed acyclic graphical models for ordered variables with non-local priors. \emph{Article submitted to Biometric Methodology}.

Shojaie, A. and Michailidis, G. (2010). Penalized likelihood methods for estimation of sparse
high-dimensional directed acyclic graphs. \emph{Biometrika} 97, 519-538.
}

\keyword{datasets}
