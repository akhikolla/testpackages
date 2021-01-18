\name{slrt}
\alias{slrt}
\title{LRT-based distances between matrices
}
\description{
  LRT-distance that we use to evaluate the performance of our covariance
  estimates. 
}
\usage{
slrt(S, trueS)
}
\arguments{
  \item{S}{estimated covariance matrix}
  \item{trueS}{true covariance matrix.}
}
\details{
    Note that this is not actually a distance in a sense that
    slrt(M1,M2) != slrt(M2,M1)}
\value{
  scalar LRT-distance
}
\references{
    Seber, G.A. (2004) Multivariate observations, Wiley
	
    Danilov, M. (2010). Robust Estimation of Multivariate Scatter under Non-Affine Equivariant Scenarios. Ph.D. thesis, 
	Department of Statistics, University of British Columbia.
}
\author{
	Mike Danilov
}
