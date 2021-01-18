\name{show.ppgasp}
\alias{show.ppgasp}
\alias{show.ppgasp-class}
\alias{show,ppgasp-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Show parllel partial Gaussian stochastic process (PP GaSP) object
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Function to print the PP GaSP model after the PP GaSP model has been constructed.
}
\usage{
\S4method{show}{ppgasp}(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ an object of  class \code{ppgasp}.}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%Provide here some details.
%}

%\references{
%% ~put references to the literature/web site here ~
%}

\author{
%%  ~~who you are~~
\packageAuthor{RobustGaSP}

Maintainer: \packageMaintainer{RobustGaSP}
}

\examples{

  library(RobustGaSP)
  
  ###PP GaSP model for the humanity model
  data(humanity_model)
  ##pp gasp
  m.ppgasp=ppgasp(design=humanity.X,response=humanity.Y,nugget.est= TRUE)
  show(m.ppgasp)
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
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
