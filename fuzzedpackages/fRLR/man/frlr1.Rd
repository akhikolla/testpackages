\name{frlr1}
\alias{frlr1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit Repeated Linear Regressions with One Variable
}
\description{
Fit a set of linear regressions which differ only in one variable.
}
\usage{
frlr1(R_X, R_Y, R_COV)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{R_X}{
the observation matrix
}
  \item{R_Y}{
the response
}
  \item{R_COV}{
common variables
}
}

\value{
the fitting results for each regression.
}
\references{
\url{https://stats.hohoweiya.xyz//regression/2017/09/26/An-R-Package-Fit-Repeated-Linear-Regressions/}
}
\author{
Lijun Wang
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
library(fRLR)
set.seed(123)
X = matrix(rnorm(50), 10, 5)
Y = rnorm(10)
COV = matrix(rnorm(40), 10, 4)
frlr1(X, Y, COV)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ repeated linear regressions }% use one of  RShowDoc("KEYWORDS")
