\name{qHardRoqF}
\alias{qHardRoqF}

\title{Hardin and Rocke F-quantiles}

\description{p-quantiles of the Hardin and Rocke (2005) scaled F distribution for squared Mahalanobis  distances based on raw MCD covariance estimators}

\usage{qHardRoqF(p, nobs, nvar, h=floor((nobs+nvar+1)/2), adj=TRUE, 
  lower.tail=TRUE, log.p=FALSE)
}

\arguments{
    \item{p}{Vector of probabilities.}
    \item{nobs}{Number of observations used in the computation of the raw MCD Mahalanobis squared distances.}
    \item{nvar}{Number of variables used in the computation of the raw MCD Mahalanobis squared distances.}
    \item{h}{Number of observations kept in the computation of the raw MCD estimate.}
    \item{adj}{logical; if TRUE (default) returns the quantile of the adjusted distribution. Otherwise returns the quantile of the asymptotic distribution.}
    \item{lower.tail}{logical; if TRUE (default), probabilities are P(X <= x) otherwise, P(X > x)}
    \item{log.p}{logical; if TRUE, probabilities p are given as log(p).}
}    

\value{The quantile of the appropriate scaled F distribution.} 

\references{
Hardin, J. and Rocke, A. (2005), The Distribution of Robust Distances.
\emph{Journal of Computational and Graphical Statistics} \bold{14}, 910--927.\cr
}

\seealso{ \code{\link{fasttle}}, \code{\link{fulltle} }


}
