\name{1-alpha quantile}
\alias{criticalvalue}


\title{The 1-alpha quantile of globaltest}

\description{Robbins and Pitman Algorithm to calculate the criticalvalue given eigenvalue vector and alpha level. }

\usage{
  criticalvalue(lam, alpha = 0.05)
}

\arguments{
  \item{lam}{The numeric vector with eigenvalues as elements.}
  \item{alpha}{The type I error rate allowed. The default is 0.05.}
}

\value{Returns a real number.}

\author{
Ningning Xu

Maintainer: Ningning Xu <n.xu@lumc.nl; xu15263142750@gmail.com>
}
\references{
   Ningning Xu, Aldo solari, Jelle Goeman, Clsoed testing with global test, with applications on metabolomics data, arXiv:2001.01541, https://arxiv.org/abs/2001.01541
}

