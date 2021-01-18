\name{MergeEdgelists}
\alias{MergeEdgelists}
\title{Combine edgelists into a single data.frame}
\usage{
  MergeEdgelists(...)
}
\arguments{
  \item{...}{data.frames, edgelists to be merged.}
}
\value{
  data.frame, single multinetwork edgelist
}
\description{
  Given two or more edgelists, create a single edgelist
  with multiple columns, two for the from and to nodes and
  one for the weights from each constituent network.
}
\examples{
edgelist1 <- data.frame(expand.grid(letters[1:2], letters[1:2]),
                        uniform=runif(4))
edgelist2 <- data.frame(v1=c("a", "a"), v2=c("a", "b"), manual=c(.3, .5))
MergeEdgelists(edgelist1, edgelist2)
}
\author{
  Stephen R. Haptonstahl \email{srh@haptonstahl.org}
}
\references{
  \url{http://www.haptonstahl.org/R}
}
\seealso{
  \code{\link{EdgelistFill}}
}

