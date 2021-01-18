\name{xtractR_M}
\alias{xtractR_M}
\title{xtractR_M}
\description{
     Internal function. Formats the output for the DetMCD algorithm.
}
\usage{
  xtractR_M(out2,X)
}
\arguments{
	\item{out2}{A list. Typically the result of a call to DetMCD_RW.}
	\item{X}{a numeric matrix or data frame without missing values.}
}
\value{
  returns a list.
}
%\examples{}
\author{
  Vakili Kaveh
}
\seealso{
  \code{\link{DetMCD}},\code{\link{DetMCD_RW}}.
}
\keyword{multivariate}
\keyword{robust}
\keyword{deterministic}
