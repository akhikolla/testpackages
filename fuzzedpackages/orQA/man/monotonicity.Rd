\name{monotonicity}
\alias{monotonicity}
\title{
Summarize results from \code{\link{pttest}} 
}
\description{
Translate results from \code{\link{pttest}} into the four categories:
overall upward trend (up), overall downward trend (down), no significant
trend (none), significant non-monotonicity (anti).
}
\usage{
monotonicity(o,alpha)
}
\arguments{
  \item{o}{
    result object from \code{\link{pttest}}
  }
  \item{alpha}{
    desired alpha level for the test
  }
}
\details{
  Translate results from \code{\link{pttest}} into the four categories:
overall upward trend (up), overall downward trend (down), no significant
trend (none), significant non-monotonicity (anti).
}
\value{
  A vector of length corresponding to the number of tests coding each
  result as either "up","down","none", or "anti".
}
\author{
  Florian Klinglmueller <float_at_lefant.net>
}
\examples{
groups <- rep(1:4,each=18)
ind <- rep(rep(1:3,each=6),4)
tdir <- sample((1:3)-2,100,rep=TRUE)
data <- matrix(rnorm(7200),nc=72)+(tdir \%*\% t(groups))
out <- pttest(data,groups,B=1000,rep=ind)
sigdir <- monotonicity(out)
table(sigdir)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }

