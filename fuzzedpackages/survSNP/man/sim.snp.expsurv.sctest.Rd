\name{sim.snp.expsurv.sctest}
\alias{sim.snp.expsurv.sctest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation of Cox Score Statistic
}
\description{
  This function simulates the  asymptotic \emph{P}-value for the Cox
  score statistic under a specified genetic risk model
}
\usage{
sim.snp.expsurv.sctest(n, gtprev, lam, a, b, ztest, diag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
    Sample size
  }
  \item{gtprev}{
    Relative genotypic frequency
  }
  \item{lam}{
    Exponential hazard rates for conditional time to event survival
    functions
  }
  \item{a}{
    Lower bound for the uniform censoring bound
  }
  \item{b}{
    Upper bound for the uniform censoring bound
  }
  \item{ztest}{
    Assumed genetic model
  }
  \item{diag}{
    Set to TRUE if print out of diagnostics is desired
  }
}
\details{
  This function simulates the genotypes 0, 1 or 2 from relative
  genotypic frequencies. Then conditional on genotype, it simulates the
  time to event from an exponential distribution with parameter
  \eqn{lambdag}{\lambda_g}. The censoring time is drawn from a uniform law on the
  interval \eqn{[a,b]}{[a,b]}.
  It is not intended to be called directly by the user.
  To conduct power calculations, use \code{sim.snp.expsurv.power}
  or the convenience wrapper function \code{survSNP.power.table}.

}
\value{
  A vector of length 2 containing the observed event rate and the
  asymptotic \emph{P}-value for the Cox score test
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
