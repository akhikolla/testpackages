\name{gen.phiMean}
\alias{gen.phiMean}
\title{Average kinship}
\description{Returns the average kinship among pairs of specified individuals}
\usage{gen.phiMean( phiMatrix)}
\arguments{
\item{phiMatrix}{A square matrix of kinship coefficients obtained for example with gen.phi. Required.}
}
\value{returns a numeric value}
\seealso{

\code{\link{gen.genealogy}}
\code{\link{gen.phi}}
\code{\link{gen.phiOver}}
\code{\link{gen.phiCI}}
}
\examples{
data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
kinship<-gen.phi(genJi)
gen.phiMean(kinship)

\donttest{
data(genea140) 
gen140<-gen.genealogy(genea140) 
phi6subjects<-gen.phi(gen140, pro=c(454422, 676521, 677273, 717537, 717634, 717709, 868572))
gen.phiMean(phi6subjects)
}
}
\keyword{manip}
