\name{surv.exp.gt.model}
\alias{surv.exp.gt.model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Exponential Hazard Rates for a Genetic Risk Model
}

\description{
  This function calculates the exponential hazard rates for
  the conditional survival functions
}
\usage{
surv.exp.gt.model(pilm, lm, gtprev, GRR, zmodel, interval)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pilm}{
    Probability that the time-to-event is greater than lm
  }
  \item{lm}{
    Landmark time used for powering the study
  }
  \item{gtprev}{
    Relative genotypic frequency
  }
  \item{GRR}{
    Genotype Hazard Ratio
  }
  \item{zmodel}{
    Genetic Risk Model (choices are "additive", "recessive" or
    "dominant")
  }
  \item{interval}{
    Interval to search for baseline hazard rate
  }
}
\details{
  Determines the time to event distributions in the three component
  mixture model discussed in the reference below.
  It is not intended to be called directly by the user.
  To conduct power calculations, use \link{sim.snp.expsurv.power}
  or the convenience wrapper function \link{survSNP.power.table}.
}
\value{
A vector of length three containing the exponential hazard rates
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
%\examples{
%# surv.exp.gt.model(pilm, lm, gtprev, GRR, model, interval)
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
