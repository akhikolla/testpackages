\name{process.data}
\alias{process.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	PROCESS DATA
}
\description{
	Initializes components necessary to perform agglomerative analysis.
}
\usage{
process.data(member, X, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{member}{Segment membership vector for the time series.}
  \item{X}{A matrix containing the time series with observations in R^d.}
  \item{alpha}{Index used for determining the distance between and
  	within segments.}
}
\details{
Called by the e.agglo method, and should not be called by the user.
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
