\name{phenoData}
\alias{phenoData}
\docType{data}
\title{Simulated Phenotype Data}
\description{
  Simulated phenotype data.
}
\usage{data(phenoData)}

\format{
  A matrix with 50 rows representing individuals and 10 phenotypes arranged in the columns.
}

\details{
  Just simulated, numerical phenotype data for testing purposes. The data consists out of 50 samples. The first 30 rows represent cases, the last 20 the controls.
  The phenotype data follows simply a N(0,1) distribution, however, the first four phenotypes were shifted by +1. 
}

\source{
  Own simulation, code can be given upon request.
}

\examples{
data(phenoData)
}
\keyword{datasets}
 
