\name{DetMCD_RW}
\alias{DetMCD_RW}
\title{DetMCD_RW}
\description{
     Internal function. Carries the re-weighting part of the DetMCD algorithm.
}
\usage{
  DetMCD_RW(ll,hlst,Xw,out2,scale_est,alpha)
}
\arguments{
	\item{ll}{integer in \code{1:6}.}
	\item{hlst}{a vector of integers between in (n/2,n).}
\item{Xw}{a n by p data matrix.}
	\item{out2}{a list. Typically the result of a call to \code{"DetMCD_CS"}.}
   \item{scale_est}{a character string specifying the
  variance functional.  Possible values are \code{"qn"}
              for the Qn or  \code{"tau"} for the tau scale.}
  \item{alpha}{a vector of values in [1/2,1].}
}
\value{
  returns a list.
}
%\examples{}
\author{
  Vakili Kaveh
}
\seealso{
  \code{\link{DetMCD}},\code{\link{DetMCD_CS}}.
}
\keyword{multivariate}
\keyword{robust}
\keyword{deterministic}
