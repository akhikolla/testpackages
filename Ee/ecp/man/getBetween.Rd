\name{getBetween}
\alias{getBetween}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GET BETWEEN DISTANCE
}
\description{
Returns the energy distance between two sets of numerical data.
}
\usage{
getBetween(alpha_, X_, Y_)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha_}{
	A weighting parameter used for calculating the energy distance. This 
	value should be in (0,2].
}
  \item{X_}{
	A n by d matrix of the n d-dimensional observations.
}
  \item{Y_}{
	A m by d matrix of the m d-dimensional observations.
}
}
\details{
The matrices X_ and Y_ do not need to have the same number of rows, but 
they do require the same number of columns.
}
\value{
The returned value is a real number indicating the energy distance between 
the two data sets.
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
\code{\link{e.agglo}}
\code{\link{e.divisive}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
set.seed(100)
X = matrix(rnorm(100),ncol=2)
Y = matrix(rnorm(126,1),ncol=2)
alpha = 1
between.distance = getBetween(alpha,X,Y)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
