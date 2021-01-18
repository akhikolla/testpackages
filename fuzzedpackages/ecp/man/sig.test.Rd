\name{sig.test}
\alias{sig.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	SIGNIFICANCE TEST
}
\description{
	Performs a permutation test.
}
\usage{
sig.test(D, R, changes, min.size, obs, env=emptyenv())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{A n by n distnace matrix.}
  \item{R}{The number of permutations to use in the permutation test.}
  \item{changes}{The set of current change points.}
  \item{min.size}{Minimum number of observations between change points.}
  \item{obs}{Test statistic value for non-permuted data.}
  \item{env}{Environment with information used to reduce computational time.}
}
\details{
Called by the e.divisive method, and should not be called by the user.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The returned value is the approximate p-value obtained by the permutation 
test. The permutaiton test is performed using the method outlined in 
Gandy (2009).
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

Gandy, A. (2009) "Sequential implementation of Monte Carlo tests with uniformly bounded resampling risk." Journal of the American Statistical Association.

Rizzo M.L., Szekely G.L. (2005). Hierarchical clustering via joint between-within distances: Extending ward's minimum variance method. Journal of Classification. pp. 151 - 183.

Rizzo M.L., Szekely G.L. (2010). Disco analysis: A nonparametric extension of analysis of variance. The Annals of Applied Statistics. pp. 1034 - 1055.
}
\author{
Nicholas A. James
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{e.divisive}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
%%\examples{
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.
%%
%%## The function is currently defined as
%%function (alpha_, X_) 
%%.Call("getWithin", alpha_, X_, PACKAGE = "energyChangePoint")
%%}% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
