\name{print.asso}
\alias{print.asso}
\title{Prints association study results}
\usage{
  \method{print}{asso} (x, ...)
}
\arguments{
  \item{x}{The association study results obtained from the
  \code{\link{AssoTestProc}}.}

  \item{...}{Usual arguments passed to the print function.}
}
\description{
  Prints formatted results from the association study
  returned by \code{\link{AssoTestProc}}.
}
\examples{
# Fit the data under the assumption that there are 3 clusters
asso.fit <- AssoTestProc(signal=signal,fam=fam,envirX=envirX,phi=phi,N=3,varSelection='PC.9')
print(asso.fit)
}
\author{
  Meiling Liu
}

