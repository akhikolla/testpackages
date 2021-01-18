\name{sigdirPttest}
\alias{sigdirPttest}
\title{
Summarize results from \code{\link{pttest}} 
}
\description{
Translate results from \code{\link{pttest}} into trinomial coding where
1 signifies an upward trend, 0 a none significant result, and -1 a
downward trend.
}
\usage{
sigdirPttest(o,alpha)
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
  Translate results from \code{\link{e2test}} into trinomial coding where
1 signifies an upward trend, 0 a none significant result, and -1 a
downward trend.
}
\value{
  a matrix with rows corresponding to the number of tests and columns
  corresponding to the number treatment groups coding each
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
out <- pttest(data,groups,B=1000,rep=ind)
sigdir <- sigdirPttest(out)
head(sigdir)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }

