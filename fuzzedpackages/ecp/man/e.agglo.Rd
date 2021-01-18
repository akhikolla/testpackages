\name{e.agglo}
\alias{e.agglo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	ENERGY AGGLOMERATIVE
}
\description{
	An agglomerative hierarchical estimation algorithm for multiple change point analysis.
}
\usage{
e.agglo(X, member=1:nrow(X), alpha=1, penalty=function(cps){0})
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{A T x d matrix containing the length T time series with d-dimensional observations.}
  \item{member}{Initial membership vector for the time series.}
  \item{alpha}{Moment index used for determining the distance between and
  	within clusters.}
  \item{penalty}{Function used to penalize the obtained goodness-of-fit statistics. This function 
  takes as its input a vector of change point locations (cps).}
}
\details{
Homogeneous clusters are created based on the initial clustering provided by the \emph{member} 
argument. In each iteration, clusters are merged so as to maximize a goodness-of-fit 
statistic. The computational complexity of this method is \emph{O(T^2)}, where \emph{T} is the 
number of observations.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a list with the following components.
\item{merged}{A (T-1) x 2 matrix indicating which segments were merged at each 
	step of the agglomerative procedure.}
\item{fit}{Vector showing the progression of the penalized goodness-of-fit 
	statistic.}
\item{progression}{A T x (T+1) matrix showing the progression of the set of change points.}
\item{cluster}{The estimated cluster membership vector.}
\item{estimates}{The location of the estimated change points.}
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
\code{\link{e.divisive}}

%% ~related works ~
Rizzo M.L., Szekely G.L. (2005). Hierarchical clustering via joint between-within distances: Extending ward's minimum variance method. Journal of Classification. pp. 151 - 183.

Rizzo M.L., Szekely G.L. (2010). Disco analysis: A nonparametric extension of analysis of variance. The Annals of Applied Statistics. pp. 1034 - 1055.
}
\examples{
set.seed(100)
mem = rep(c(1,2,3,4),times=c(10,10,10,10))
x = as.matrix(c(rnorm(10,0,1),rnorm(20,2,1),rnorm(10,-1,1)))
y = e.agglo(X=x,member=mem,alpha=1,penalty=function(cp,Xts) 0)
y$estimates


\dontrun{
# Multivariate spatio-temporal example
# You will need the following packages:
#	mvtnorm, combinat, and MASS
library(mvtnorm); library(combinat); library(MASS)
set.seed(2013)
lambda = 1500 #overall arrival rate per unit time
muA = c(-7,-7) ; muB = c(0,0) ; muC = c(5.5,0)
covA = 25*diag(2)
covB = matrix(c(9,0,0,1),2)
covC = matrix(c(9,.9,.9,9),2)
time.interval = matrix(c(0,1,3,4.5,1,3,4.5,7),4,2)
#mixing coefficents
mixing.coef = rbind(c(1/3,1/3,1/3),c(.2,.5,.3), c(.35,.3,.35), 
	c(.2,.3,.5))
stppData = NULL
for(i in 1:4){
	count = rpois(1, lambda* diff(time.interval[i,]))
	Z = rmultz2(n = count, p = mixing.coef[i,])
	S = rbind(rmvnorm(Z[1],muA,covA), rmvnorm(Z[2],muB,covB),
		rmvnorm(Z[3],muC,covC))
	X = cbind(rep(i,count), runif(n = count, time.interval[i,1],
		time.interval[i,2]), S)
	stppData = rbind(stppData, X[order(X[,2]),])
}
member = as.numeric(cut(stppData[,2], breaks = seq(0,7,by=1/12)))
output = e.agglo(X=stppData[,3:4],member=member,alpha=1,
	penalty=function(cp,Xts) 0)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ agglomerative }
\keyword{ hierarchical }% __ONLY ONE__ keyword per line
