\name{e.split}
\alias{e.split}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	ENERGY SPLIT
}
\description{
	Finds the most likely location for a change point across all current 
	clusters.
}
\usage{
	e.split(changes, D, min.size, for.sim=FALSE, env=emptyenv())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{changes}{
		A vector containing the current set of change points.
	}
	\item{D}{
		An n by n distance matrix.
	}
	\item{min.size}{
		Minimum number of observations between change points.
	}
	\item{for.sim}{
		Boolean value indicating if the function is to be run on permuted 
		data for significance testing.
	}
	\item{env}{
		Environment that contains information to help reduce computational
		time.
	}
}
\details{
This method is called by the e.divisive method, and should not be called by 
the user.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	A list with the following components is returned.
	\item{first}{The index of the first element of the cluster to be divided.}
	\item{second}{The index of the last element of the cluster to be divided.}
	\item{third}{The new set of change points.}
	\item{fourth}{The distance between the clusters created by the newly 
		proposed change point.
	}
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
