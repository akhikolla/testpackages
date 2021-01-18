\name{plot.mdr}
\docType{methods}
\alias{plot.mdr}
\alias{plot,mdr-method}

\title{Plot an MDR Object}
\description{
 The function offers informative plots for an \code{mdr} object.
}
\usage{
 \method{plot}{mdr}(x, which=NULL, ...)
}

\arguments{
  \item{x}{Object of class \code{mdr}.}
  \item{which}{Specifies for which association fold should the plot be created.}
  \item{...}{Additional plotting parameters.}
}

\details{
 This function plots the density of the precision of an \code{mdr} object.
}

\author{Daniel Fischer}

\examples{
# Perform the MDR
# res <- mdr(X=temp,status=status,fold=3,top=20)

# plot(res)
}

\keyword{methods}
\keyword{hplot}
