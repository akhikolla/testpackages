\name{kcpa}
\alias{kcpa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Kernel Change Point Analysis
}
\description{
	An algorithm for multiple change point analysis that uses the 'kernel trick' and dynamic programming.
}
\usage{
kcpa(X, L, C)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
  	A T x d matrix containing the length T time series with d-dimensional observations.
  	}
  \item{L}{
	The maximum number of change points.
	}
  \item{C}{
	The constant used to penalize the inclusion of additional change points in the fitted model.
	}
}
\details{
Segments are found through the use of dynamic programming and the kernel trick.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	If the algorithm determines that the best fit is obtained through using k change points then the returned value is an array of length k, containing the change point locations.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
Arlot S., Celisse A., Harchaoui Z. (2019). A Kernel Multiple Change-point Algorithm via Model Selection. J. Mach. Learn. Res., 20, 162:1-162:56.
}
\author{
Nicholas A. James
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\keyword{kernel}% __ONLY ONE__ keyword per line

