\name{quanff}
\alias{quanff}
\title{quanff}
\description{
     Internal function. Converts alpha values to h values.
}
\usage{
  quanff(alpha,n,p)
}
\arguments{
  \item{alpha}{a value in [1/2,1].}
  \item{n,p}{integers.}
}
\value{
  returns an integer.
}
%\details{}
\references{
   Hubert, M., Rousseeuw, P.J. and Verdonck, T. (2012),
   "A deterministic algorithm for robust location and scatter", Journal of
   Computational and Graphical Statistics, in press.
}
\examples{
quanff(0.75,n=100,p=5);
}
\author{
  Vakili Kaveh
}
\keyword{multivariate}
\keyword{robust}
\keyword{deterministic}
