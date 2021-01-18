\name{inQn}
\alias{inQn}
\title{Test function for the qn}
\usage{
  inQn(x)
}
\arguments{
 \item{x}{Vector of 2 or more numbers. Should contain no ties.}
}
\value{
  the value of the qn estimator of scale.
}
\description{
  Test function for the qn used in DetR.
}
\references{
    see \code{pcaPP::qn} and citation("pcaPP").
}
\examples{
set.seed(123)	#for reproductibility
x<-rnorm(101)
inQn(x)
#should be the same:
pcaPP::qn(x)
}
\author{
Kaveh Vakili. Calls code translated from the cde for computing the Qn found in 
 package pcaPP (by  Peter Filzmoser, Heinrich Fritz, Klaudius Kalcher , see citation("pcaPP")).
}
\keyword{multivariate}
\keyword{robust}

