\name{print.eqtl}
\docType{methods}
\alias{print.eqtl}
\alias{print-method}
\alias{print,eqtl-method}

\title{Print an eqtl Object}
\description{
  Prints an \code{eqtl} object.
}
\usage{
 \method{print}{eqtl}(x, which=NULL, sig=0.01, output="bed", ...)
}

\arguments{
  \item{x}{Object of class \code{eqtl}.}
  \item{which}{Which center gene should be printed.}
  \item{sig}{Significance level.}
  \item{output}{Output format.}
  \item{...}{Additional parameters.}
}

\details{
  The function prints SNPs in the surroundings of a gene from an \code{eqtl} object. 

  By default all genes are considered, subsets can be defined with the \code{which}
  option. The \code{sig} option gives the threshold which results should be shown.
}

\author{Daniel Fischer}

\examples{
\dontrun{
myeqtl <- eQTL(geneMatrix,genoData,singleLoc,genoSamples,singleSamples,windowSize,method="LM")
myeqtl
print(myeqtl, sig=0.05)
}
}

\keyword{methods}
\keyword{print}
