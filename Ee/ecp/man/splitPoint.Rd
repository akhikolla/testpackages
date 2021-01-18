\name{splitPoint}
\alias{splitPoint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	SPLIT POINT
}
\description{
	Finds the most likely location of a change point within a given 
	segment.
}
\usage{
splitPoint(start, end, D, min.size)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{start}{the index of the first observation in a segment.}
  \item{end}{The index of the last observation in a segment.}
  \item{D}{A n by n distance matrix.}
  \item{min.size}{Minimum number of observations between change points.}
}
\details{
Called by the e.divisive method, and should not be called by the user.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The returned value is a vector. The first component is the most likely 
position of a change point. The second component is the distance between 
the segments created by this proposed change point.
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
%%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
