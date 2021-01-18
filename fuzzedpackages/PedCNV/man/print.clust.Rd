\name{print.clust}
\alias{print.clust}
\title{Prints clustering results}
\usage{
  \method{print}{clust} (x, ...)
}
\arguments{
  \item{x}{The clustering results obtained from the
  \code{\link{ClusProc}}.}

  \item{...}{Usual arguments passed to the print function.}
}
\description{
  Prints formatted results returned by
  \code{\link{ClusProc}}.
}
\examples{
# Fit the data under the given clustering numbers
clus.fit <- ClusProc(signal=signal,N=2:6,varSelection='PC.9')
print(clus.fit)
}
\author{
  Meiling Liu and Sungho Won
}

