\name{environ.4.data}
\alias{environ.4.data}
\alias{environ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Environmental model function 
%%  ~~function to do ... ~~
}
\description{

This function is the environmental model where the output is the concentration of the pollutant at the space-time grid used in in Bliznyuk et al. (2008). The physical inputs are 4 dimensions.
%It is an example of a nonpolynomial model which exhibits a shape similar to that of a multivariate polynomial. 

}
\usage{
environ.4.data(x, s=c(0.5, 1, 1.5, 2, 2.5), t=seq(from=0.3, to=60, by=0.3))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a 4-dimensional vector specifying input parameters.}
%%     ~~Describe \code{x} here~~
  \item{s}{  spatial location to be evaluated.}
  %%     ~~Describe \code{s} here~~
  \item{t}{  time point to be evaluated.}
  %%     ~~Describe \code{s} here~~

}

\details{
Bliznyuk, N., Ruppert, D., Shoemaker, C., Regis, R., Wild, S., & Mugunthan, P. (2008). Bayesian calibration and uncertainty analysis for computationally expensive models using optimization and radial basis function approximation. \emph{Journal of Computational and Graphical Statistics}, \bold{17}(2), 270-294.

}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A real value of this function evaluated at \code{x}.

}
\references{
%% ~put references to the literature/web site here ~

N. Bliznyuk, D. Ruppert, C. Shoemaker, R. Regis, S. Wild, & P. Mugunthan (2008). Bayesian calibration and uncertainty analysis for computationally expensive models using optimization and radial basis function approximation. \emph{Journal of Computational and Graphical Statistics}, \bold{17}(2), 270-294.


S. Surjanovic, D. Bingham, Virtual Library of Simulation Experiments:  Test Functions
and Datasets, retrieved March 29, 2016, from \url{http://www.sfu.ca/~ssurjano/fried.html}.


}

\author{
\packageAuthor{RobustGaSP}

Maintainer: \packageMaintainer{RobustGaSP}
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

%\examples{
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.
%%
%%## The function is currently defined as
%%function (x) 
%%{
%%  }
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{optimize}
\keyword{internal}% __ONLY ONE__ keyword per line
