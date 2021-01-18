\name{e.divisive}
\alias{e.divisive}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	ENERGY DIVISIVE
}
\description{
	A divisive hierarchical estimation algorithm for multiple change point analysis.
}
\usage{
e.divisive(X, sig.lvl=.05, R=199, k=NULL, min.size=30, alpha=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
  	A T x d matrix containing the length T time series with d-dimensional observations.
  	}
	\item{sig.lvl}{
		The level at which to sequentially test if a proposed change point is 
		statistically significant.
	}
	\item{R}{
		The maximum number of random permutations to use in each iteration of the permutation test. The 
		permutation test p-value is calculated using the method outlined in Gandy (2009).
	}
	\item{k}{
		Number of change point locations to estimate, suppressing permutation based testing. 
		If k=NULL then only the statistically significant estimated change points are returned.
	}
	\item{min.size}{
		Minimum number of observations between change points.
	}
	\item{alpha}{
		The moment index used for determining the distance between and within 
		segments.
	}
}
\details{
Segments are found through the use of a binary bisection method and a permutation 
test. The computational complexity of this method is \emph{O(kT^2)}, where \emph{k} is the 
number of estimated change points, and \emph{T} is the number of observations.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	The returned value is a list with the following components.
	\item{k.hat}{The number of clusters within the data created 
		by the change points.}
	\item{order.found}{The order in which the change points were estimated.}
	\item{estimates}{Locations of the statistically significant change points.}
	\item{considered.last}{Location of the last change point, that 
		was not found to be statistically significant at the given significance level.}
	\item{permutations}{The number of permutations performed by each of the sequential permutation test.}
	\item{cluster}{The estimated cluster membership vector.}
	\item{p.values}{Approximate p-values estimated from each permutation test.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
Matteson D.S., James N.A. (2013). A Nonparametric Approach for Multiple Change Point Analysis of Multivariate Data.

Nicholas A. James, David S. Matteson (2014). "ecp: An R Package for Nonparametric 
Multiple Change Point Analysis of Multivariate Data.", "Journal of Statistical Software, 
62(7), 1-25", URL "http://www.jstatsoft.org/v62/i07/"
}
\author{
Nicholas A. James
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{e.agglo}}

%% ~related works ~
Gandy, A. (2009) "Sequential implementation of Monte Carlo tests with uniformly bounded resampling risk." Journal of the American Statistical Association.

Rizzo M.L., Szekely G.L (2005). Hierarchical clustering via joint between-within distances: Extending ward's minimum variance method. Journal of Classification.

Rizzo M.L., Szekely G.L. (2010). Disco analysis: A nonparametric extension of analysis of variance. The Annals of Applied Statistics.
}
\examples{
set.seed(100)
x1 = matrix(c(rnorm(100),rnorm(100,3),rnorm(100,0,2)))
y1 = e.divisive(X=x1,sig.lvl=0.05,R=199,k=NULL,min.size=30,alpha=1)
x2 = rbind(MASS::mvrnorm(100,c(0,0),diag(2)),
	MASS::mvrnorm(100,c(2,2),diag(2)))
y2 = e.divisive(X=x2,sig.lvl=0.05,R=499,k=NULL,min.size=30,alpha=1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ divisive }
\keyword{ hierarchical }% __ONLY ONE__ keyword per line
