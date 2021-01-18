\encoding{UTF-8}
\name{dataset}
\docType{data}
\alias{dataset}
\title{
 simulated data
}
\description{
  This data file consists of six simulated predictors or variables with three class categories. For each class category the values are independently generated from the normal distribution with the mean µ1, µ2 and µ3 and the variances held at unity. The means are varied such that the problems range from near-separable problems, to near-random.
}
\usage{dataset}

\format{A data.frame containing 300 observations of six variables.}

\source{
  Landgrebe T, Duin R (2006) A simplified extension of the Area under the ROC to the multiclass domain. In: Proceedings 17th Annual Symposium of the Pattern Recognition Association of South Africa. PRASA, pp. 241–245.
}

\seealso{
 \code{\link{sim}}
}

\examples{
# load the dataset
data(dataset)
}

\keyword{datasets}
