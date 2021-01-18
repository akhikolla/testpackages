\encoding{UTF-8}
\name{data_test}
\docType{data}
\alias{data_test}
\title{
 simulated data
}
\description{
  This data file consists of 300 objects with 10 features. The features x1-x5 are informative and define the cluster structure of the dataset. The clusters are generated in the two-dimensional space x1-x2. The values of the features x3-x5 are identically generated as for x2. Features values x1-x5 are normally distributed values with the same standard deviation and different mean values. Features x6-x10 are random variables uniformly distributed in the interval [0, 1] and present the uninformative features.
}
\usage{data_test}

\format{A data.frame containing 300 observations of 11 variables and class with three labels.}

\source{
Landgrebe T, Duin R (2006) A simplified extension of the Area under the ROC to the multiclass domain. In: Proceedings 17th Annual Symposium of the Pattern Recognition Association of South Africa. PRASA, pp. 241–245.
}

\seealso{
  \code{\link{leukemia72}}, \code{\link{datasetF6}}, \code{\link{leukemia72_2}}, \code{\link{leukemia_miss}}
}

\examples{
# load the dataset
data(data_test)
}

\keyword{datasets}
