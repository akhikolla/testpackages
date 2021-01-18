\name{multidcov}
\alias{multidcov}
\title{
  Symmetric multivariate distance covariance}
\description{
  Calculate either the symmetric or asymmetric multivariate
  distance covariance statistic.
}
\usage{
multidcov(S,symmetric=TRUE,alpha=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{the n x d matrix for which you wish to calculate the 
  dependence between d columns from n samples}
  \item{alpha}{A scaling parameter in the interval (0,2] used for calculating distances.}
  \item{symmetric}{logical; if TRUE (the default), calculates the symmetric version of the multivariate distance covariance. See details.}
  }
\details{
If symmetric==TRUE, calculates: sum_{i=1}^{d} dcovustat(S[,i],S[,-i])
If symmetric==FALSE, calculates: sum_{i=1}^{d-1} dcovustat(S[,i],S[,(i+1):d])
}
\value{
returns a scalar equal to the multivariate distance covariance statistic
for the columns of S
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
David Matteson
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
nObs <- 1024
nComp <- 3

simM <- matrix(rnorm(nComp*nComp),nComp)

# simulate some data:
simS<-cbind(rgamma(nObs, shape = 1, scale = 2),
            rgamma(nObs, shape = 3, scale = 2),
            rgamma(nObs, shape = 9, scale = 0.5))


simS <- scale(simS) #Standardize variance for identifiability

#mix the sources:
xData <- simS \%*\% simM

multidcov(simS) #close to zero
multidcov(whitener(xData)$Z) #should be larger than simS
multidcov(xData) #greater than zero
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distance }
\keyword{ covariance }% __ONLY ONE__ keyword per line
\keyword{ dcov}
