\name{matern_5_2_deriv}
\alias{matern_5_2_deriv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
The derivative of matern correlation function with regard to the inverse range parameter }
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function computes the derivative of a correlation matrix using 1-dimensional matern correlation function, with regard to the inverse range parameter.
}
\usage{
matern_5_2_deriv(R0_i, R, beta_i)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{R0_i}{
%%     ~~Describe \code{R0_i} here~~
an absolute difference matrix of the i-th input vector.

}
  \item{R}{
%%     ~~Describe \code{R} here~~
the correlation matrix. 
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
A matrix in which each element is the derivative of 1-dimensional matern correlation function with regard to the inverse range parameter.
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

\seealso{
\code{\link{matern_5_2_funct}}.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{internal}
