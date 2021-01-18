\name{sim.snp.expsurv.power}
\alias{sim.snp.expsurv.power}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Asymptotic and Empirical Power
}
\description{
This function calculates asymptotic and empirical power for SNP
association studies
}
\usage{
sim.snp.expsurv.power(GHR, B, n, raf, erate, pilm, lm, model, test, alpha, 
                      exactvar = FALSE, interval = c(0, 10), rootint = c(0.1, 200))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GHR}{
    Genotype Hazard Ratio
  }
  \item{B}{
    Number of simulation replicates (set to 0 if no empirical
    calculations are desired)
  }
  \item{n}{
    Sample size
  }
  \item{raf}{
    Relative risk allele frequency
  }
  \item{erate}{
    Event Rate
  }
  \item{pilm}{
    Probability that the time-to-event is greater than lm
  }
  \item{lm}{
    Landmark time used for powering the study
  }
  \item{model}{
    True genetic risk model (choices are "additive", "recessive" or
    "dominant". For the asymptotic calculations only the "additive"
    should be used until further notice)
  }
  \item{test}{
    Hypothesized genetic risk model (choices are "additive", "recessive"
    or "dominant". For the asymptotic calculations only the "additive" should
    be used until further notice)
  }
  \item{alpha}{
    Nominal two-sided type I error rate
  }
  \item{exactvar}{
    Indicator for using the exact variance formula
  }
  \item{interval}{
    Interval to search for baseline hazard rate
  }
  \item{rootint}{
    Interval to search for censoring bound
  }
}
\details{
This function calculates asymptotic and empirical power for SNP
association studies
}
\value{
  This function returns a \code{data.frame} with the following columns
  \item{B}{
    Number of simulation replicates
  }
  \item{raf}{
    Relative allelic frequency of the B allele
  }
  \item{q0}{
    Relative genotypic frequency for AA
  }
  \item{q1}{
    Relative genotypic frequency for AB
  }
  \item{q2}{
    Relative genotypic frequency for BB
  }
  \item{lam0}{
    Exponential hazard rate for P(T>t|AA)
  }
  \item{lam1}{
    Exponential hazard rate for P(T>t|AB)
  }
  \item{lam2}{
    Exponential hazard rate for P(T>t|BB)
  }
  \item{GHR}{
    Genotype Hazard Ratio
  }
  \item{pilm}{
    Probability that the time to event in the population
    exceeds the landmark lm: P(T>lm)
  }
  \item{lm}{
    Landmark (see pilm above)
  }
  \item{model}{
    The true genetic risk model
  }
  \item{test}{
    The hypothesized genetic risk model
  }
  \item{a}{
    Parameter for the censoring distribution (uniform on [a,b])
  }
  \item{b}{
    Parameter for the censoring distribution (uniform on [a,b])
  }
  \item{erate}{
    Event rate
  }
  \item{n}{
    Sample size
  }
  \item{powB}{
    Empirical Power
  }
  \item{pow}{
    Asymptotic Power based on the exact variance formula
  }
  \item{pow0}{
    Asymptotic Power based on the approximate variance formula
  }
  \item{v1}{
    First term of the variance (v1+v2+v12)
  }
  \item{v2}{
    Second term of the variance (v1+v2+v12)
  }
  \item{v12}{
    Third term (covariance) of the variance (v1+v2+v12)
  }
  
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

\seealso{
\code{\link{survSNP.power.table}}
}
\examples{
# See vignette for examples including details on the example
# considered below

results<-sim.snp.expsurv.power(GHR=1.25, B=0, n=500, raf=0.1, erate=0.75, pilm=0.5, 
                               lm=1, model="additive", test="additive", alpha=0.05)
results[,c("n","erate","alpha","pow0")]
}

