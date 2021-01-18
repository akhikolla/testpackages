\name{tt.edges}
\alias{tt.edges}
\docType{data}
\title{
   Pathogen link over the time
}
\description{
  A data.frame of all transmission links between hosts (pathogens).
  Five columns ID, source, target, time and weight.
}
\usage{data("transmissiontree")}
\format{
  A data frame with 13 observations on the following 5 variables.
  \describe{
    \item{\code{ID}}{a numeric vector}
    \item{\code{source}}{a character vector}
    \item{\code{target}}{a factor with levels \code{113} \code{104} \code{116} \code{115} \code{111} \code{109} \code{105} \code{108} \code{106} \code{112}}
    \item{\code{time}}{a character vector}
    \item{\code{weight}}{a character vector}
  }
}
\examples{
data(transmissiontree)
print(tt.edges)
}
\keyword{datasets}
