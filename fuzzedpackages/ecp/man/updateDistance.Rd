\name{updateDistance}
\alias{updateDistance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	UPDATE DISTANCE
}
\description{
	Update the distance between the newly created segment and the other 
	segments.
}
\usage{
updateDistance(i,j,K,ret)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{i}{The segment that makes up the left portion of the new segment.}
  \item{j}{The segment that makes up the right portion of the new segment.}
  \item{K}{Integer indicating the progress of the agglomerative process.}
  \item{ret}{A list with 'gof', 'list', 'N', 'sizes', 'right', 'left', 
  	'open', 'D', and 'lm' components.}
}
\details{
This method is called by the e.agglomerative method, and should not be called 
by the user.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a list with the following components.
\item{gof}{Vector showing the progression of the goodness of fit statistic.}
\item{list}{Matrix showing the progression of the set of change points.}
\item{N}{Number of initial segments.}
\item{sizes}{Sizes of each segment during the agglomerative process.}
\item{right}{Vector containing indices of the right adjacent segments.}
\item{left}{Vector containing indices of the left adjacent segments.}
\item{open}{Vector indicating if a segment has been merged.}
\item{D}{Matrix of distances between segments.}
\item{lm}{Vector containing indices of the starting point of a segment.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
James NA, Matteson DS (2013). A Nonparametric Approach for Multiple Change Point Analysis of Multivariate Data.

Nicholas A. James, David S. Matteson (2014). "ecp: An R Package for Nonparametric 
Multiple Change Point Analysis of Multivariate Data.", "Journal of Statistical Software, 
62(7), 1-25", URL "http://www.jstatsoft.org/v62/i07/"

Rizzo ML, Szekely GL (2005). Hierarchical clustering via joint between-within distances: Extending ward's minimum variance method. Journal of Classification. pp. 151 - 183.
}
\author{
Nicholas A. James
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{e.agglo}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
