\name{gmultidcov}
\alias{gmultidcov}
\title{
  Symmetric multivariate distance covariance for grouped components}
\description{
  Calculate either the symmetric or asymmetric multivariate
  distance covariance statistic for a given grouping of the components.
}
\usage{
gmultidcov(S,group=1:ncol(S),alpha=1,symmetric=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{The n x d matrix for which you wish to calculate the 
  dependence between d columns from n samples}
  \item{group}{A length d vector which indicates group membership for each component}
  \item{alpha}{A scaling parameter in the interval (0,2] used for calculating distances.}
  \item{symmetric}{logical; if TRUE (the default), calculates the symmetric version of the multivariate distance covariance. See details.}
  }
\details{
Suppose that the groups are numbered 1,2,...,C and that group is a vector indicating 
group membership for each component.
If symmetric==TRUE, calculates: sum_{i=1}^{C} dcovustat(S[,group==i],S[,group!=i])
If symmetric==FALSE, calculates: sum_{i=1}^{C-1} dcovustat(S[,group==i],S[,group>i])
}
\value{
Returns a scalar equal to the multivariate distance covariance statistic for grouped 
components of S.
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

\examples{

library(steadyICA)
S = matrix(rnorm(300),ncol=3)
group = c(1,2,2)
gmultidcov(S,group,TRUE) # close to zero
gmultidcov(S,group,FALSE) # sill close to zero

Sigma = matrix(c(1,0.7,0,0.7,1,-0.2,0,-0.2,1),ncol=3)
X = MASS::mvrnorm(100,rep(0,3),Sigma)
gmultidcov(X,group,TRUE) # further from zero
gmultidcov(X,group,FALSE) # further from zero

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distance }
\keyword{ covariance }% __ONLY ONE__ keyword per line
\keyword{ dcov}
