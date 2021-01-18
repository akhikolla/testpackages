\name{gof.update}
\alias{gof.update}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	GOODNESS OF FIT UPDATE
}
\description{
	Updates the goodness of fit statistic.
}
\usage{
gof.update(i,ret)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{i}{segments which is to be merged with the segment that is adjacent 
  	to its right.}
  \item{ret}{A list with 'gof', 'right', 'left', 'D', and 'size' comonents.}
}
\details{
Called by the e.agglo method, and should not be called by the user.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a real number. This is the updated goodness of fit statistic.
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
