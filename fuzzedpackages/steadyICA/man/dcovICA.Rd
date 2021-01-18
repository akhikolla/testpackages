\name{dcovICA}
\alias{dcovICA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ICA via distance covariance for 2 components}
\description{
This algorithm finds the rotation which minimizes the distance covariance between two orthogonal components via the angular parameterization of a 2x2 orthogonal matrix with the function stats::optimize. The results will be (approximately)
equivalent to steadyICA but this function is much faster (but does not extend to higher dimensions).}
\usage{
dcovICA(Z, theta.0 = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{The whitened n x d data matrix, where n is the number of observations and d the number of components.}
  \item{theta.0}{Determines the interval to be searched by the optimizer: lower bound = theta.0, upper
  bound = pi/2. Changing theta.0 affects the initial value,
  where the initial value = theta.0+(1/2+sqrt(5)/2)*pi/2, see \link{optimize}.}
  }
  
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{theta.hat}{Estimated minimum.}
\item{W}{W = t(theta2W(theta.hat))} 
\item{S}{Estimated independent components.}
\item{obj}{The distance covariance of S.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Matteson, D. S. & Tsay, R. Independent component analysis via U-Statistics. 
    <http://www.stat.cornell.edu/~matteson/#ICA>
}

\author{
David Matteson and Benjamin Risk
}

\seealso{
\code{\link{steadyICA}},
\code{\link{optimize}}
}
\examples{
library(JADE)
library(ProDenICA)
set.seed(123)
simS = cbind(rjordan(letter='j',n=1024),rjordan(letter='m',n=1024))
simM = mixmat(p=2)
xData = simS\%*\%simM
xWhitened = whitener(xData)

#Define true unmixing matrix as true M multiplied by the estimated whitener:
#Call this the target matrix:
W.true <- solve(simM\%*\%xWhitened$whitener) 


a=Sys.time()
est.dCovICA = dcovICA(Z = xWhitened$Z,theta.0=0)
Sys.time()-a

#See the example with steadyICA for an explanation
#of the parameterization used in amari.error:
amari.error(t(est.dCovICA$W),W.true)

##NOTE: also try theta.0 = pi/4 since there may be local minima
  \dontrun{est.dcovICA = dcovICA(Z = xWhitened$Z,theta.0=pi/4)
  amari.error(t(est.dcovICA$W),W.true)}

a=Sys.time()
est.steadyICA = steadyICA(X=xWhitened$Z,verbose=TRUE)
Sys.time()-a
amari.error(t(est.steadyICA$W),W.true)
##theta parameterization with optimize is much faster
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ givens }
\keyword{ distance }% __ONLY ONE__ keyword per line
\keyword{ covariance }
\keyword{ ICA }