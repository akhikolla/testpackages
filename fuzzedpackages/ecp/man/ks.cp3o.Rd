\name{ks.cp3o}
\alias{ks.cp3o}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	CHANGE POINTS ESTIMATION BY PRUNED OBJECTIVE (VIA KOLMOGOROV-SMIRNOV STATISTIC)
}
\description{
	An algorithm for multiple change point analysis that uses dynamic programming and 
	pruning. The Kolmogorov-Smirnov statistic is used as the goodness-of-fit measure.
}
\usage{
ks.cp3o(Z, K=1, minsize=30, verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{Z}{
  		A T x d matrix containing the length T time series with d-dimensional observations.
  	}
	\item{K}{
		The maximum number of change points.
	}
	\item{minsize}{
		The minimum segment size.
	}
	\item{verbose}{
		A flag indicating if status updates should be printed.
	}
}
\details{
Segmentations are found through the use of dynamic programming and 
pruning. For long time series, consider using ks.cp3o_delta.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	The returned value is a list with the following components.
	\item{number}{The estimated number of change points.}
	\item{estimates}{The location of the change points estimated by the procedure.}
	\item{gofM}{A vector of goodness of fit values for differing number of 
		change points. The first entry corresponds to when there is only 
		a single change point, the second for when there are two, and 
		so on.}
	\item{cpLoc}{The list of locations of change points estimated by the procedure for 
		different numbers of change points up to K.}		
	\item{time}{The total amount to time take to estimate the change point locations.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
W. Zhang, N. A. James and D. S. Matteson, "Pruning and Nonparametric Multiple Change Point Detection," 2017 IEEE International Conference on Data Mining Workshops (ICDMW), New Orleans, LA, 2017, pp. 288-295.
}
\author{
Wenyu Zhang
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%% ~related works ~
Kifer D., Ben-David S., Gehrke J. (2004). Detecting change in data streams. International Conference on Very Large Data Bases.
}
\examples{
set.seed(400)
x = matrix(c(rnorm(100),rnorm(100,3),rnorm(100,0,2)))
y = ks.cp3o(Z=x, K=7, minsize=30, verbose=FALSE)
#View estimated change point locations
y$estimates
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ pruning }
\keyword{ dynamic }% __ONLY ONE__ keyword per line
