\name{asypow}
\alias{asypow}

\title{
Calculating the asymptotic power and variance
}
\description{
  This function calculates the asymptotic power and variance
  assuming that the survival distribution is a mixture of
  exponentials with rates and the censoring distribution is
  uniform on the interval (a,b).
}
\usage{
asypow(n, theta, a, b, lambda0, q, p, alpha, z,exactvar)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
    Sample size
  }
  \item{theta}{
    Effect size (log genotype hazard ratio (GHR))
  }
  \item{a}{
    Censoring distribution parameter (assumed to be uniform on [a,b])
  }
  \item{b}{
    Censoring distribution parameter (assumed to be uniform on [a,b])
  }
  \item{lambda0}{
    Baseline exponential hazard rate
  }
  \item{q}{
    Relative risk allele frequency
  }
  \item{p}{
    Relative genotype frequency 
  }
  \item{alpha}{
    Nominal two-sided type I error rate
  }
  \item{z}{
    Genotype scores (right now only additive scores AA=0,AB=1,BB=2
    generate correct power)
  }
  \item{exactvar}{
    Indicator for using the exact variance formula
  }
}
\details{ This function is called by \link{sim.snp.expsurv.power}
  to calculate the asymptotic variance (exact and approximate)
  formulas. It is not intended to be called directly by the user.
  To conduct power calculations, use \link{sim.snp.expsurv.power}
  or the convenience wrapper function \link{survSNP.power.table}.
}
\value{

  \item{power}{Asymptotic power based on exact variance formula}
  \item{power0}{Asymptotic power based on approximate variance formula}
  \item{v1}{First term of asymptotic variance}
  \item{v2}{Second term of asymptotic variance}
  \item{v12}{Third term of the asymptotic variance (covariance)}
  \item{vapprox}{Approximate asymptotic variance formula (=v1)}
  \item{exact}{Exact asymptotic variance formula (=v1+v2+v12)}
  \item{diff}{Difference between variances (=v2+v12)}
  \item{ratio}{Ratio of variances (=v1/(v1+v2+v12))}
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
