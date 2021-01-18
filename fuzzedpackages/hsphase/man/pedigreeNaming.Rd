\name{pedigreeNaming}
\alias{pedigreeNaming}
\title{
Fix Pedigree Errors
}
\description{
Tries to link the inferred pedigree from \code{\link{rpoh}} with the sire IDs in the original pedigree and fix pedigree errors.}
\usage{
pedigreeNaming(inferredPedigree, realPedigree)
}
\arguments{
  \item{inferredPedigree}{
{\code{matrix} inferred pedigree  (output of \code{\link{rpoh}} )}

}

\item{realPedigree}{
{\code{matrix} } original pedigree

}

}

\details{
This function calls the \code{\link{bmh}} and \code{\link{recombinations}} functions to count the number of recombinations in each half-sib group.  
}
\value{
Returns the inferred pedigree with the best fit to the sire names used in the original pedigree file.
}

\seealso{
 \code{\link{rpoh}} and \code{\link{ohg}} 
}
\examples{
# Please run demo(hsphase)
}

\keyword{pedigree}
\keyword{infer}


