\name{splitPointC}
\alias{splitPointC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
SPLIT POINT-C
}
\description{
C++ function that is called by splitPoint() to perform calculations.
}
\usage{
splitPointC(s_, e_, D_, min_size_)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s_}{
	Index of the first obervation in a segment.
}
  \item{e_}{
	Index of the last observation in a segment.
}
  \item{D_}{
	A distance matrix.
}
  \item{min_size_}{
	The minimum segment size.
}
}
\details{
As with the splitPoint method, this method should not be calle by the user.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	Returns a vector. The first component is the most likely position of 
	a change point. The second component is the distance between the 
	segments created by this proposed change point.

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

Rizzo ML, Szekely GL (2010). Disco analysis: A nonparametric extension of analysis of variance. The Annals of Applied Statistics. pp. 1034 - 1055.
}
\author{
Nicholas A. James
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{splitPoint}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
