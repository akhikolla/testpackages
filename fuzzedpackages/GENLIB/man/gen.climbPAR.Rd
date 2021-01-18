\name{gen.climbPAR}
\alias{gen.climbPAR}
\title{Internal function used by gen.findMRCA}
\description{gen.climbPAR is an internal function used by gen.findMRCA.}
\usage{gen.climbPAR(gen, individuals, founder)}
\arguments{
\item{gen}{The genealogy to consider.}
\item{individuals}{probands to consider}
\item{founder}{the founder by which the distance is calculated}
}
\value{returns a list containing the founder ID and the distance.}
\seealso{
\code{\link{gen.findMRCA}}
\code{\link{gen.getAncestorsPAR}}
}
\examples{
data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
gen.findMRCA(genJi, individuals=c(1,29), NbProcess = 1)
}
\keyword{methods}
