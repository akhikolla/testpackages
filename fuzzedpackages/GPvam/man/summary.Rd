\name{summary}
\alias{summary.GPvam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary
}
\description{
Prints summary information for object of class GPvam
}
\usage{
\method{summary}{GPvam}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
object of class GPvam
}
  \item{\dots}{
other arguments to be passed to summary
}
}


\author{
Andrew Karl \email{akarl@asu.edu}
Yan Yang
Sharon Lohr
}


\seealso{
 \code{\link{plot.GPvam}}
}
\examples{
\dontrun{
data(vam_data)
result<-GPvam(vam_data)
summary(result)
}
}
\keyword{ regression }