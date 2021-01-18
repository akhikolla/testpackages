\name{frlr2}
\alias{frlr2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit Repeated Linear Regressions with Two Variables
}
\description{
Fit a set of linear regressions which differ only in two variables.
}
\usage{
frlr2(R_X, R_idx1, R_idx2, R_Y, R_COV)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{R_X}{
the observation matrix
}
  \item{R_idx1}{
the first identical feature
}
  \item{R_idx2}{
the second identical feature
}
  \item{R_Y}{
the response variable
}
  \item{R_COV}{
common variables
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
the fitting results for each regression.
}
\references{
\url{https://stats.hohoweiya.xyz//regression/2017/09/26/An-R-Package-Fit-Repeated-Linear-Regressions/}
}
\author{
Lijun Wang
}

\examples{
library(fRLR)
set.seed(123)
X = matrix(rnorm(50), 10, 5)
Y = rnorm(10)
COV = matrix(rnorm(40), 10, 4)
idx1 = c(1, 2, 3, 4, 1, 1, 1, 2, 2, 3)
idx2 = c(2, 3, 4, 5, 3, 4, 5, 4, 5, 5)
frlr2(t(X), idx1, idx2, Y, t(COV))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ repeated linear regressions }% use one of  RShowDoc("KEYWORDS")
