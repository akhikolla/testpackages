\name{censbnd}
\alias{censbnd}

\title{
Calculate bound for censoring distribution.
}
\description{
This function computes the bound 
}
\usage{
censbnd(lambda, p, crate, rootint = c(0.1, 1000))
}
\arguments{
  \item{lambda}{
    Baseline exponential hazard rate
  }
  \item{p}{
    Relative genotype frequency
  }
  \item{crate}{
    Desired censoring rate
  }
  \item{rootint}{
    Interval to be searched for the root
  }
}
\details{
  The time to event distribution \eqn{T}{T} is assumed to be a mixture of
  exponentials with parameter \eqn{lambda0}{\lambda_0}, \eqn{lambda1}{\lambda_1}
  and \eqn{lambda2}{\lambda_2} with mixing proportion
  \eqn{p0}{p_0},\eqn{p1}{p_1} and \eqn{p2}{p_2}.
  Suppose that the censoring distribution is uniform on the interval
  \eqn{[0,b]}{[0,b]}. 
  This function calculates \eqn{b}{b} for a desired censoring rate.
  It is not intended to be called directly by the user.
  To conduct power calculations, use \link{sim.snp.expsurv.power}
  or the convenience wrapper function \link{survSNP.power.table}.
}
\value{
This function returns a list from the \code{uniroot} function. The
root component of this list in the  
}
\references{
  Kouros Owzar, Zhiguo Li, Nancy Cox and Sin-Ho Jung.
  Power and Sample Size Calculations for SNP Association
  Studies with Censored Time-to-Event Outcomes.
  \url{http://onlinelibrary.wiley.com/doi/10.1002/gepi.21645/abstract}
}
\author{
Kouros Owzar, Zhiguo Li, Nancy Cox, Sin-Ho Jung and Chanhee Yi
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
censbnd(0.1,hwe(0.1),0.9)$root
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
