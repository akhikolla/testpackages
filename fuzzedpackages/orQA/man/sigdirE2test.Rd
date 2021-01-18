\name{sigdirE2test}
\alias{sigdirE2test}
\title{
Summarize results from \code{\link{e2test}} 
}
\description{
Translate results from \code{\link{e2test}} into trinomial coding where
1 signifies an upward trend, 0 a none significant result, and -1 a
downward trend.
}
\usage{
sigdirE2test(o,alpha)
}
\arguments{
  \item{o}{
    result object from \code{\link{e2test}}
  }
  \item{alpha}{
    desired alpha level for the test
  }
}
\details{
  Translate results from \code{\link{e2test}} into trinomial coding where
1 signifies an upward trend, 0 a none significant result, and -1 a
downward trend.
}
\value{
  a vector of length corresponding to the number of tests coding each
  result as either 1,0, or -1.
}
\author{
  Florian Klinglmueller <float_at_lefant.net>
}
\examples{
groups <- rep(1:4,each=18)
ind <- rep(rep(1:3,each=6),4)
tdir <- sample((1:3)-2,100,rep=TRUE)
data <- matrix(rnorm(7200),nc=72)+(tdir \%*\% t(groups))
out <- e2test(data,groups,B=1000,rep=ind)
sigdir <- sigdirE2test(out)
table(sigdir)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }

