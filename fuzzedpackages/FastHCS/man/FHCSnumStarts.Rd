\name{FHCSnumStarts}
\alias{FHCSnumStarts}
\title{Computes the number of starting q-subsets}
\usage{
  FHCSnumStarts(q,gamma=0.99,eps=0.5)
}
\arguments{
  \item{q}{Number of desired components for the PCA model.}

  \item{gamma}{Desired probability of having at least one 
	clean starting q-subset.}

  \item{eps}{suspected contamination rate of the sample.}

}
\value{
  An integer number of starting q-subsets.
}
\description{
  Computes the number of starting q-subsets to take so that there is a 99% chance that at least one is outlier free.
  This is an internal function not intended to be called 
  by the user.
}
\examples{
FHCSnumStarts(q=3,gamma=0.99,eps=0.4)
}
\author{
Kaveh Vakili
}
\keyword{multivariate}
\keyword{robust}
