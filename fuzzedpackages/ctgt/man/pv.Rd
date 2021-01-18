\name{p-value}
\alias{pv}

\title{p-value of globaltest}

\description{Robbins and Pitman Algorithm to calculate the p-value given the observed value and the eigenvalue vector.}

\usage{
  pv(x, lam)
}

\arguments{
  \item{x}{The observed value that is used to calculate the coresponding the p-value.}
  \item{lam}{The numeric vector with eigenvalues as elements.}
}

\value{Returns a value between 0 and 1. }

\author{
Ningning Xu

Maintainer: Ningning Xu <n.xu@lumc.nl; xu15263142750@gmail.com>
}
\references{
   Ningning Xu, Aldo solari, Jelle Goeman, Clsoed testing with global test, with applications on metabolomics data, arXiv:2001.01541, https://arxiv.org/abs/2001.01541
}

