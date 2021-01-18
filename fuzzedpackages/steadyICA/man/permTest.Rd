\name{permTest}
\alias{permTest}
\title{
  Permutation test for mutual independence.
}
\description{
  Calculates an approximate p-values based upon a permutation test for mutual 
  independence.
}
\usage{
permTest(S, group=1:ncol(S), R=199, FUN=c('gmultidcov','compInd'), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{The n x d matrix for which you wish to test the 
  dependence between d columns from n samples}
  \item{group}{A length d vector which indicates group membership for each component}
  \item{R}{The number of permutations to perform in order to obtain the approximate p-value.}
  \item{FUN}{The function used to determine mutual independence. This is one of either gmultidcov 
  or compInd.}
  \item{...}{Additionl arguments passed to FUN. See details.}
  }
\details{
Suppose that the groups are numbered 1,2,...,C and that group is a vector indicating 
group membership for each component.
If symmetric==TRUE, calculates: sum_{i=1}^{C} dcovustat(S[,group==i],S[,group!=i])
If symmetric==FALSE, calculates: sum_{i=1}^{C-1} dcovustat(S[,group==i],S[,group>i])

If no additional arguments are supplied for FUN then the default values are used. In the case of gmultidcov, values for alpha and symmetric can be supplied. While for compInd only the value of alpha is needed.
}
\value{
Returns an approximate p-values based upon a permutation test.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Nicholas James
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{dcovustat}},
\code{\link[energy:dcov]{energy::dcov}}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\examples{
%
%
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distance }
\keyword{ covariance }% __ONLY ONE__ keyword per line
\keyword{ dcov}
