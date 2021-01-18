\name{gen.getAncestorsPAR}
\alias{gen.getAncestorsPAR}
\title{Internal function used by gen.findMRCA}
\description{gen.getAncestorsPAR is an internal function used by gen.findMRCA.}
\usage{gen.getAncestorsPAR(gen, pro)}
\arguments{
\item{gen}{The genealogy to consider.}
\item{pro}{probands to consider}
}
\value{Vector of the founders identity.}
\seealso{
\code{\link{gen.findMRCA}}
\code{\link{gen.climbPAR}}
}
\examples{
data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
gen.findMRCA(genJi, individuals=c(1,29), NbProcess = 1)
}
\keyword{methods}
