\name{gen.getFoundersPAR}
\alias{gen.getFoundersPAR}
\title{Internal function used by gen.findFounders}
\description{gen.getFoundersPAR is an internal function used by gen.findFounders.}
\usage{gen.getFoundersPAR(gen, pro)}
\arguments{
\item{gen}{The genealogy to consider.}
\item{pro}{probands to consider}
}
\value{A vector of the founders identity.}
\seealso{
\code{\link{gen.findFounders}}
\code{\link{gen.findDistance}}
}
\examples{
data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
gen.founder(genJi, individuals=c(1,29))
gen.findFounders(genJi, individuals=c(1,29), NbProcess = 1)
}
\keyword{methods}
