\name{simM3PL}
\alias{simM3PL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Simulate M3PL item responses
}
\description{
  Simulates item responses according to the Multidimensional 3-Parameter Logistic (M3PL) model
}
\usage{
  simM3PL(ipar, cors, p, n.simulee = 100, D = 1, easiness = T, seed = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ipar}{a data frame containing M3PL item parameters, specifically a1, a2, ... , d, and c}
  \item{cors}{a square matrix of the lower diagonal elements of a variance-covariance (VCV) matrix, 
   including 1's in the main diagonal}
  \item{p}{the number of latent dimensions}
  \item{n.simulee}{the number of simulees to generate}
  \item{D}{scaling constant: 1.7 or 1.0 (default=1.0)}
  \item{easiness}{logical, T if d is related to the \emph{easiness} of items per Reckase, F otherwise}
  \item{seed}{random number seed}
}
\details{
  This function simulates item responses according to the Multidimensional 3-Parameter Logistic (M3PL) model using
  the item parameters input to the function. Thetas are drawn from the multivariate standard normal distribution with 
  the population variance-covariance (correlation) matrix input to the function.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{call}{function call stack}
  \item{theta}{a \emph{n.simulee} by \emph{p} matrix of true theta values}
  \item{resp}{a data frame of simulated item responses named "R1", "R2", ...}
%% ...
}
\references{
  Reckase, M. D. (2009). Multidimensional Item Response Theory. New York: Springer.
}
\author{
  Seung W. Choi
}

\examples{
  data(sample.ipar)
  vcv1<-diag(3)
  vcv1[lower.tri(vcv1,diag=FALSE)]<-c(.5,.6,.7)
  resp1<-simM3PL(sample.ipar, vcv1, 3, n.simulee = 100, seed = 1234)$resp
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Psychometrics }
\keyword{ IRT }

