\name{dcovustat}
\alias{dcovustat}
\title{
  Calculate distance covariance via U-statistics}
\description{
  Calculates the square of the U-statistic formulation of distance covariance. This is  faster than the function 'dcov' in the R package 'energy'  and requires less memory. Note that negative values are possible in this version.}
\usage{
dcovustat(x,y,alpha=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A vector or matrix.}
  \item{y}{A vector or matrix with the same number of observations as x, though the number of columns of x and y may differ}
  \item{alpha}{A scaling parameter in the interval (0,2] used for calculating distances.}
  }
\value{
Returns the distance covariance U-statistic.
}
\references{
Matteson, D. S. & Tsay, R. Independent component analysis via U-Statistics. 
    <http://www.stat.cornell.edu/~matteson/#ICA>

Szekely, G., Rizzo, M. & Bakirov, N. Measuring and testing dependence by correlation of distances. (2007) \emph{The Annals of Statistics}, 35, 2769-2794.
}
\author{
David Matteson
}
\note{
The value returned by dcovustat is equal to the square of the value returned by energy::dcov in the limit. 

In dcovustat, a vector of length n is stored; in energy::dcov, an n x n matrix is stored. Thus, dcovustat requires far less memory and works for very large datasets.

Even though dcovustat converges to the square of the distance covariance of the random variables x and y, it can be negative.
}

\seealso{
\code{\link{multidcov}},
\code{\link[energy:dcov]{energy::dcov}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
x = rnorm(5000)
y = rbinom(5000,1,0.5)
y = y - 1*(y==0)
z = y*exp(-x) #some non-linear dependence

dcovustat(x[1:1000],y[1:1000]) #close to zero

a = Sys.time()
dcovustat(x[1:1000],z[1:1000]) #greater than zero
a = Sys.time() - a

#measures of linear dependence close to zero:
cov(x,z)
cor(rank(x),rank(z))


\dontrun{
#dcovustat differs from energy::dcov but are equal in the limit
library(energy)
b = Sys.time()
(dcov(x[1:1000],z[1:1000]))^2
b = Sys.time() - b
as.double(b)/as.double(a) #dcovustat is much faster

## energy::dcov and dcovustat become approximately equal as n increases:
c = Sys.time()
dcovustat(x,z)
c = difftime(Sys.time(), c, sec)
d = Sys.time()
(dcov(x,z)^2)
d = difftime(Sys.time(), d, sec)
as.double(d)/as.double(c) 
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distance }
\keyword{ covariance }% __ONLY ONE__ keyword per line
\keyword{ dcov}
