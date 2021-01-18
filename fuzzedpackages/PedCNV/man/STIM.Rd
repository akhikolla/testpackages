\name{STIM}
\alias{STIM}
\title{Score test with the intensity value}
\usage{
  STIM(envirX, signal, fam, alpha, phi, sig2g, sig2)
}
\arguments{
  \item{envirX}{The matrix of environmental variables. The
  intercept should be included if it's needed.}

  \item{fam}{The FAM file which follows the format defined
  in PLINK.}

  \item{signal}{The matrix of intensity measurements. The
  row names must be consistent with the Individual ID in
  fam file.}

  \item{alpha}{The estimated parameters for environmental
  variables under null hypothesis. This value can be
  calculated by using function \code{\link{AssoTestProc}}.}

  \item{phi}{The matrix of correlation between
  individuals.}

  \item{sig2g}{The estimated standard error for polygenic
  effect under null hypothesis. This value can be
  calculated by using function \code{\link{AssoTestProc}}.}

  \item{sig2}{The estimated standard error for
  environmental effect under null hypothesis. This value
  can be calculated by using function
  \code{\link{AssoTestProc}}.}
}
\value{
  It returns the statistic value and pvalue of the score
  test. \item{STIMs}{The statistic value of score test with
  the intensity value under null hypothesis.}
  \item{STIMp}{The pvalue of score test with the intensity
  value under null hypothesis.} \item{df}{The degree of
  freedom of score test with the intenstiy value under null
  hypothesis.}
}
\description{
  Calculates the score test statistics with the intensity
  value.
}
\examples{
# Fit the data under the assumption that there are 3 clusters
asso.fit <- AssoTestProc(signal=signal,fam=fam,envirX=envirX,phi=phi,N=3,varSelection='PC.9')
alpha <- asso.fit$para$alpha
sig2g <- asso.fit$para$sig2g
sig2 <- asso.fit$para$sig2
STIM(envirX=envirX,signal=signal,fam=fam,alpha=alpha,phi=phi,sig2g=sig2g,sig2=sig2)
}
\author{
  Meiling Liu, Sungho Won
}

