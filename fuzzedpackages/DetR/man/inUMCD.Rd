\name{inUMCD}
\alias{inUMCD}
\title{Test function for unimcd}
\usage{
  inUMCD(x)
}
\arguments{
 \item{x}{Vector of 2 or more numbers. Should contain no ties.}
}
\value{
  the value of the unimcd estimator of scale.
}
\description{
  Test function for the unimcd used in DetR.
}
\references{
     Rousseeuw, P. J. (1984), Least Median of Squares Regression, Journal of the American Statistical Association,79,
 871--880. 
}
\examples{
set.seed(123)	#for reproductibility
x<-rnorm(101)
inUMCD(x)
}
\author{
Kaveh Vakili
}
\keyword{multivariate}
\keyword{robust}

