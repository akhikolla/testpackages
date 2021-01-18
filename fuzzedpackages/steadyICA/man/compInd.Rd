\name{compInd}
\alias{compInd}
\title{
  Complete Measure of Mutual Multivariate Independence}
\description{
  Calculates a complete empirical measure of mutual multivariate independence. Makes 
  use of the utils::combn function.
}
\usage{
compInd(S,group=1:ncol(S),alpha=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{The n x d matrix for which you wish to calculate the 
  dependence between d columns from n samples.}
  \item{group}{A length d vector which indicates group membership for each component.}
  \item{alpha}{The index used in calculating the distance between sample observations.}
  }
%\details{
%}
\value{
Returns a scalar equal to the empirical multivariate distance between the observed samples, 
and their grouped counterpart.}
\references{
%% ~put references to the literature/web site here ~
Chasalow, Scott (2012) combinat: Combinatorics Utilities <http://CRAN.R-project.org/package=combinat
}
\author{
Nicholas James
}
\note{
Suppose that the each component belongs to exactly one of C groups. 
This method makes use of the utils::combn and combinat::permn functions. As a result it will be both computationally and memory intensive, even for small to moderate n and small C.
}

\seealso{
\code{\link{dcovustat}},
\code{\link[energy:dcov]{energy::dcov}}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

library(steadyICA)
library(combinat)
set.seed(100)
S = matrix(rnorm(40),ncol=4)
group = c(1,2,3,3)
compInd(S,group,1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distance }
\keyword{ covariance }% __ONLY ONE__ keyword per line
\keyword{ dcov }
