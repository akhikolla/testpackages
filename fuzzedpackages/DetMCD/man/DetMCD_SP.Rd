\name{DetMCD_SP}
\alias{DetMCD_SP}
\title{DetMCD_SP}
\description{
     Internal function. Computes the starting points for the DetMCD algorithm.
}
\usage{
  DetMCD_SP(Data,scale_est,tol)
}
\arguments{
	\item{Data}{a numeric matrix or data frame without missing values.}
   \item{scale_est}{a character string specifying the
  variance functional.  Possible values are \code{"qn"}
              for the Qn or  \code{"tau"} for the tau scale.}
   \item{tol}{a small positive numeric value to be
  used for determining numerical 0.}
}
\value{
  returns a list.
}
%\examples{}
\author{
  Vakili Kaveh
}
\seealso{
  \code{\link{DetMCD}}.
}
\keyword{multivariate}
\keyword{robust}
\keyword{deterministic}
