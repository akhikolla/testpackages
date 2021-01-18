\name{GetSampleFromFile}
\alias{GetSampleFromFile}
\title{Sample from the rows of a (possibly large) text file (NOT IMPLEMENTED)}
\usage{
  GetSampleFromFile(n, out.of, filename)
}
\arguments{
  \item{n}{numeric, size of sample to be taken.}

  \item{out.of}{numeric, number of rows in the data set not
  including the header.}

  \item{filename}{character, name of the file containing
  the data. This must be a tab-delimited file with a header
  row formatted per the default options for
  \code{\link{read.delim}}.}
}
\value{
  data.frame, size n random subset of the rows of filename
}
\description{
  Read a large text file in batches, keeping the rows to be
  included in the sample. Return a data.frame whose rows
  are the sample.
}
\examples{
\dontrun{x <- dils:::GetSampleFromFile(10, 150, "folder/containing/data.txt")}
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
  \code{\link{GetSampleFromDb}}
}

