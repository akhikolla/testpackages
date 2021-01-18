\name{estim_herit}
\alias{estim_herit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of heritability in linear mixed models.
}
\description{
The function computes heritability when the random effects follow either a Gaussian distribution or a mixture of a Dirac and a Gaussian distribution.}
\usage{
estim_herit(Y, Z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
vector of observations of size n
}
  \item{Z}{
matrix with n rows and N columns
}
}

\value{
\item{heritability}{Heritability}
}
\references{
The method is developped in the paper "Heritability estimation in high dimensional linear mixed models" (A.Bonnet,vE.Gassiat, C.Levy-Leduc,2014)
}
\author{
Anna Bonnet
}
