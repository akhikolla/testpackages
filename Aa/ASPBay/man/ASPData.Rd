\name{ASPData}
\alias{ASPData}
\docType{data}
\title{ Simulated dataset }
\description{
  Simulations of 1000 controls and 1000 affected sib pairs with 22 SNPs. There is one causal SNP with an OR of 2.
}
\usage{data(ASPData)}
\format{
  A list with the following elements.
  \describe{
    \item{\code{Control}}{21 genotypes of controls no including causal SNP}
    \item{\code{Index}}{21 genotypes of index cases no including causal SNP}
    \item{\code{IBD}}{Vector of IBD states for each affected sib pairs}
	\item{\code{Causal}}{The name of the causal SNP}
  }
}
\details{See \code{vignette("ASPBay")}.}
\examples{
data(ASPData)
}
\keyword{datasets}
