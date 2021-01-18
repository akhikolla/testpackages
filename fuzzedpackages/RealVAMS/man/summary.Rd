\name{summary}
\alias{summary.RealVAMS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary
}
\description{
Prints summary information for object of class RealVAMS
}
\usage{
\method{summary}{RealVAMS}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
object of class RealVAMS
}
  \item{\dots}{
other arguments to be passed to summary
}
}


\author{
Andrew Karl \email{akarl@asu.edu}, Jennifer Broatch, Jennifer Green
}


\seealso{
 \code{\link{plot.RealVAMS}}
}
\examples{
\dontrun{
res<-RealVAMS(example.score.data,example.outcome.data)
summary(res)
}
}
\keyword{ regression }