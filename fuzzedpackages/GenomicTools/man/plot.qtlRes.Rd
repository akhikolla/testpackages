\name{plot.qtlRes}
\docType{methods}
\alias{plot.qtlRes}
\alias{plot,qtl-method}

\title{Plot an qtlRes Object}
\description{
 The function offers informative plots for an \code{qtlRes} object.
}
\usage{
 \method{plot}{qtlRes}(x, which=NULL, sig=0.01, verbose=TRUE,
                       log=FALSE, genome=NULL, pch=1, ...)
 }

\arguments{
  \item{x}{Object of class \code{qtlRes}.}
  \item{which}{Specifies for which phenotypes the plot should be created.}
  \item{sig}{Chosen significance level.}
  \item{verbose}{Logical, extended feedback of the function.}
  \item{log}{Logical, y-axis scale is log(base=10)-scaled.}
  \item{genome}{Data frame with Chromosome names and lengths.}
  \item{pch}{Plot symbol to by used.}
  \item{...}{Additional plotting parameters.}
}

\details{
 This function plots the test results of an \code{qtlRes} object as a Manhattan plot. In case only an \code{qtlRes} object is provided, then the genome information
 are estimated from this object. Alternatively, a data frame with genome information can be provided. In this case \code{genone} has to be a data.frame with the
 two columns \code{chr} and \code{length}, indicating the Chromosome names (similar as the ones given in \code{qtlRes}) and the corresponding lengths. The information
 for the human genome, Ensembl build 68, is also included can be used with using the option \code{genome="Human68"}.
}

\author{Daniel Fischer}

\examples{

\dontrun{
data("phenoData")
data("genotData")

qtlresult <- QTL(pheno=phenoData, geno=genotData, which = 3:4)

plot(qtlresult, which=2)

}
}

\keyword{methods}
\keyword{hplot}
