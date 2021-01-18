\name{compPcaParams}
\alias{compPcaParams}
\title{Computes the center vector, eigenvalues and loading matrix corresponding to a PCA model of a data matrix with respect to a subset of observations in a data set}
\usage{
  compPcaParams(x,fitd,q=NULL,z0=NULL,seed=1)
}
\arguments{
  \item{x}{A data matrix x.}
  \item{fitd}{The (internal) result of a call to FastHCS.}
  \item{q}{Desired rank of the SVD decomposition.}
	\item{z0}{Optional. Result of a call to \code{FHCSkernelEVD}.}
  \item{seed}{Seed used to initialize the RNG. Defaults to 1.}
}
\value{
   A  list with the following components:
	  \item{center}{The multivariate mean of the observations with indexes in best.}
	  \item{loadings}{The (rank q) loadings matrix of the observations with indexes in best.}
	  \item{eigenvalues}{The eigenvalues of the observations with indexes in best multiplied by a consistency factor.}
	  \item{scores}{The value of the projected on the space of the principal components data (the centred  data multiplied by the loadings matrix) is returned. Hence, cov(scores) is the diagonal matrix diag(eigenvalues).}
}
\description{
  This function is used in FastHCS to compute the parameter estimates of the PCA models used at different steps of the algorithm.
  It is an internal function not intended to be called 
  by the user. 
}
\author{
Kaveh Vakili, Eric Schmitt
}
\keyword{multivariate}
\keyword{robust}

