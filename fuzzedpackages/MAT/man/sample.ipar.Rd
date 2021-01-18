\name{sample.ipar}
\alias{sample.ipar}
\docType{data}
\title{
  Sample item parameters
}
\description{
  A sample item parameter file containing 180 Multidimensional 3-PL (M3PL) model. 
}
\usage{data(sample.ipar)}
\format{
  A data frame with item parameters for 180 items.
  \describe{
    \item{\code{a1}}{the discrimination parameter for theta 1}
    \item{\code{a2}}{the discrimination parameter for theta 2}
    \item{\code{a3}}{the discrimination parameter for theta 3}
    \item{\code{d}}{the easiness parameter, d=-a*b}
    \item{\code{c}}{the guessing parameter}
  }
}
\details{
  First 60 items are primarily loaded on theta 1, second 60 on theta 2, and last 60 on theta 3.
}
\examples{
  data(sample.ipar)
}
\keyword{datasets}
