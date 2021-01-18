\name{DetMCD_CS}
\alias{DetMCD_CS}
\title{DetMCD_CS}
\description{
     Internal function. Computes the Csteps for the DetMCD algorithm.
}
\usage{
  DetMCD_CS(Data,scale_est,h,out1)
}
\arguments{
	\item{Data}{a numeric matrix or data frame without missing values.}
   \item{scale_est}{a character string specifying the
  variance functional.  Possible values are \code{"qn"}
              for the Qn or  \code{"tau"} for the tau scale.}
   \item{h}{a vector of integers (between n/2 and n).}
	\item{out1}{A list. Typically the result of a call to DetMCD_SP.}
}
\value{
  returns a list.
}
%\examples{}
\author{
  Vakili Kaveh
}
\seealso{
  \code{\link{DetMCD}},\code{\link{DetMCD_SP}}.
}
\keyword{multivariate}
\keyword{robust}
\keyword{deterministic}
