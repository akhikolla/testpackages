\name{perm.cluster}
\alias{perm.cluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	PERMUTE CLUSTERS
}
\description{
	Permutes time series observations within specified segments.
}
\usage{
perm.cluster(D, points)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{A n by n distance matrix.}
  \item{points}{The set of current change points.}
}
\details{
Called by the e.divisive method, and should not be called by the user.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns the n by n distance matrix for the permuted data.
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

Rizzo ML, Szekely GL (2010). Disco analysis: A nonparametric extension of analysis of variance. The Annals of Applied Statistics. pp. 1034 - 1055.
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
