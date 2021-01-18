\name{GetSampleFromDataFrame}
\alias{GetSampleFromDataFrame}
\title{Randomly select rows from a data.frame.}
\usage{
  GetSampleFromDataFrame(n, x)
}
\arguments{
  \item{n}{numeric, size of sample.}

  \item{x}{data.frame, data whose rows will be sampled.}
}
\value{
  data.frame, size n random subset of the rows of x
}
\description{
  Randomly select \code{n} rows from data.frame \code{x}.
}
\examples{
data(iris)   # provides example data
x <- dils:::GetSampleFromDataFrame(10, iris)
}
\author{
  Stephen R. Haptonstahl \email{srh@haptonstahl.org}
}
\references{
  \url{https://github.com/shaptonstahl/}
}
\seealso{
  \code{\link{ScalablePCA}},
  \code{\link{GetSampleFromFile}},
  \code{\link{GetSampleFromFile}}
}

