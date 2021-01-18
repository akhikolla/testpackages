\name{GetSampleFromDb}
\alias{GetSampleFromDb}
\title{Sample from the rows of a (possibly large) database table (NOT IMPLEMENTED)}
\usage{
  GetSampleFromDb(n, db)
}
\arguments{
  \item{n}{numeric, size of sample to be taken.}

  \item{db}{connection, connection to the database table
  containing the data.}
}
\value{
  data.frame, size n random subset of the rows of filename
}
\description{
  Access a database table directly. Return a data.frame
  whose rows are the sample.
}
\examples{
\dontrun{x <- dils:::GetSampleFromDb(10, my.db)}
}
\author{
  Stephen R. Haptonstahl \email{srh@haptonstahl.org}
}
\references{
  \url{https://github.com/shaptonstahl/}
}
\seealso{
  \code{\link{ScalablePCA}},
  \code{\link{GetSampleFromDataFrame}},
  \code{\link{GetSampleFromFile}}
}

