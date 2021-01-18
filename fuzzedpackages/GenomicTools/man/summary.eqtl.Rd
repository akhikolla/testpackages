\name{summary.eqtl}
\docType{methods}
\alias{summary.eqtl}
\alias{summary-method}
\alias{summary,eqtl-method}

\title{Summarize an eqtl Object}
\description{
  Summarizes and prints an \code{eqtl} object in an informative way.
}
\usage{
 \method{summary}{eqtl}(object, sig=0.01, ...)
}

\arguments{
  \item{object}{Object of class \code{eqtl}.}
  \item{sig}{Significance level to print.}
  \item{...}{Additional parameters.}
}

\details{
 This function gives a summary of an \code{eqtl} object. 
}

\author{Daniel Fischer}

\examples{
# First perform an eQTL
# lm.myEQTL <- eQTL(gex=ourExpression, geno=genotData, xAnnot=xAnnotDF, windowSize=1,
#                   genoSamples=genoSamples, method="LM")

# summary(lm.myEQTL)
}

\keyword{methods}
\keyword{print}
