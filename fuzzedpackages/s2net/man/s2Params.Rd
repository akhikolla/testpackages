\name{s2Params}
\alias{s2Params}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Hyper-parameter wrapper for \code{s2net}
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This is a very simple function that collapses the input parameters into a named vector to supply to C++ methods. 
}
\usage{
s2Params(lambda1, lambda2 = 0, gamma1 = 0, gamma2 = 0, gamma3 = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lambda1}{
%%     ~~Describe \code{lambda1} here~~
elastic-net regularization parameter - \eqn{l_1} norm.
}
  \item{lambda2}{
%%     ~~Describe \code{lambda2} here~~
elastic-net regularization parameter - \eqn{l_2} norm.

}
  \item{gamma1}{
%%     ~~Describe \code{gamma1} here~~
s2net weight hyper-parameter.
}
  \item{gamma2}{
%%     ~~Describe \code{gamma2} here~~
s2net covariance hyper-parameter (between 1 and \code{Inf}).
}
  \item{gamma3}{
%%     ~~Describe \code{gamma3} here~~
s2net shift hyper-parameter (between 0 and 1).
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns a named vector of S3 class \code{s2Params}.
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
%%\author{
%%  ~~who you are~~
%%}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{s2Data}}, \code{\link{s2Fista}}
}
% R documentation directory.
%%\keyword{  }% use one of  RShowDoc("KEYWORDS")
