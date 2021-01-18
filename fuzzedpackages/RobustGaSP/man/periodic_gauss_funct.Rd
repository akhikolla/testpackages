\name{periodic_gauss_funct}
\alias{periodic_gauss_funct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Periodic folding of the Gaussian correlation function . }
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function computes the values of periodic folding of the Gaussian correlation function (i.e. the power expoential correlation with roughness parameter being 2).
}
\usage{
periodic_gauss_funct(d, beta_i)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
%%     ~~Describe \code{d} here~~
locations the periodic Gaussian correlation function are to be evaluated.
}
  \item{beta_i}{
%%     ~~Describe \code{beta_i} here~~
the inverse range parameter.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A matrix in which each element is the value of the periodic folding of the Gaussian correlation function  evaluated at that location. 
}
%\references{
%% ~put references to the literature/web site here ~
%}
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

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{internal}