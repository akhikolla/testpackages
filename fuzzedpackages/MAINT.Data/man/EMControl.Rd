\name{EMControl}
\alias{EMControl}

\title{ Constructor function for objects of class EMControl}
\description{
  This function will create a control object of class \code{\linkS4class{EMControl}} containing the control parameters 
  for the EM algorithm  used in estimation of Gaussian mixtures by function \code{\link{Idtmclust}}.  
}

\usage{
  EMControl(nrep=0, maxiter=1000, convtol=0.01, protol=1e-3, seed=NULL, pertubfct=1, 
   k2max=1e6, MaxVarGRt=1e6)
}

\arguments{
  \item{nrep}{Number of replications (different randomly generated starting points) of the EM algorithm.}
  \item{maxiter}{Maximum number of iterations in each replication of the EM algorithm.}
  \item{convtol}{Numeric tolerance for testing the convergence of the EM algorithm. Convergence is assumed when the log-likelihood changes less than convtol.}
  \item{protol}{Numeric tolerance for the mixture proportions. Proportions below protol, considered to be zero, are not allowed.}
  \item{seed}{Starting value for random generator.}
  \item{pertubfct}{Perturbation factor used to control the degree similarity between the alternative randomly generated starting points of the EM algorithm. Increasing (decreasing) the value of pertubfct increases (decreases) the expected difference between the starting points generated.}
  \item{k2max}{Maximal allowed l2-norm condition number for correlation matrices. Solutions in which any component has correlation matrix with condition number above k2max, are considered to be spurious solutions and are eliminated from the EM search.}
  \item{MaxVarGRt}{Maximal allowed ratio of variances across components. Solutions in which any variable has a ratio between its maximal and minimal (across components) variances above MaxVarGRt, are considered to be spurious solutions and are eliminated from the EM search.}
}

\value{An \code{EMControl} object}

\seealso{
  \code{\link{Idtmclust}}
}

