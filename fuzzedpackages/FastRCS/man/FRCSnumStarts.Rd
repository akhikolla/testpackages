\name{FRCSnumStarts}
\alias{FRCSnumStarts}
\title{Computes the number of starting p-subsets}
\usage{
  FRCSnumStarts(p,gamma=0.99,eps=0.5)
}
\arguments{
  \item{p}{number of dimensions of the data matrix X.}

  \item{gamma}{desired probability of having at least one 
	clean starting p-subset.}

  \item{eps}{suspected contamination rate of the sample.}

}
\value{
  An integer number of starting p-subsets.
}
\description{
Computes the number of starting p-subsets so that the desired probability of selecting at least one clean one is achieved.
  This is an internal function not intended to be called 
  by the user.
}
\examples{
FRCSnumStarts(p=3,gamma=0.99,eps=0.4)
}
\author{
Kaveh Vakili
}
\keyword{regression}
\keyword{robust}



