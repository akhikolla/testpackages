\name{ScalablePCA}
\alias{ScalablePCA}
\title{Perform Principal Component Analysis on a large data set}
\usage{
  ScalablePCA(x, filename = NULL, db = NULL,
    subsample = 10000, n.subsamples = 1000, ignore.cols,
    use.cols, return.sds = FALSE, progress.bar = FALSE)
}
\arguments{
  \item{x}{data.frame, data over which to run PCA}

  \item{filename}{character, name of the file containing
  the data. This must be a tab-delimited file with a header
  row formatted per the default options for
  \code{\link{read.delim}}.}

  \item{db}{Object type, database connection to table
  containing the data (NOT IMPLEMENTED).}

  \item{subsample}{numeric or logical, If an integer, size
  of each subsample.  If FALSE, runs PCA on entire data
  set.}

  \item{n.subsamples}{numeric, number of subsamples.}

  \item{ignore.cols}{numeric, indices of columns not to
  include.}

  \item{use.cols}{numeric, indices of columns to use.}

  \item{return.sds}{logical, if TRUE return the standard
  deviations of each network's edge weights.}

  \item{progress.bar}{logical, if TRUE then progress in
  running subsamples will be shown.}
}
\value{
  If \code{return.sds} is FALSE, return named vector of
  component weights for first dimension of principal
  component analysis (see example for comparison to
  \code{\link{prcomp}}).

  If \code{return.sds} is TRUE, return a list.
  \tabular{ll}{ coefficients \tab named vector of the
  component weights for first dimension of principal
  component analysis (see example for comparison to
  \code{\link{prcomp}}).\cr sds \tab named vector of the
  standard deviations of each network's edge weights.\cr }
}
\description{
  Run \code{prcomp} on subsamples of the data set and
  compile the results for the first dimension.
}
\details{
  Scales the function \code{\link{prcomp}} to data sets
  with an arbitrarily large number of rows by running
  \code{prcomp} on repeated subsamples of the rows.
}
\examples{
data(iris)        # provides example data
prcomp(iris[,1:4], center=FALSE, scale.=FALSE)$rotation[,1]
ScalablePCA(iris, subsample=10, use.cols=1:4)
ScalablePCA(iris, subsample=10, ignore.cols=5)
}
\author{
  Stephen R. Haptonstahl \email{srh@haptonstahl.org}
}
\references{
  \url{https://github.com/shaptonstahl/}
}
\seealso{
  \code{\link{prcomp}}
}

