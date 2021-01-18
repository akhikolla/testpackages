\name{exportPLINKSet}
\alias{exportPLINKSet}
\docType{package}
\title{
  exportPLINKSet
}
\description{
  Simple function using Rcpp to write the gene set to a file in the PLINK set format.
}
\usage{
  exportPLINKSet(geneSets, fname)	
}
\arguments{
  \item{geneSets}{An object created by the \code{makeGeneSet()} function.
}
  \item{fname}{The name of the PLINK file to be created.
}
}
\value{
  A Boolean indicating if the file was successfully written.
}
\seealso{
    \code{\link{makeGeneSet}}
}
\examples{
    # Please see the vignette or the package description 
    # for an example of using this function.
}
\keyword{ file }
