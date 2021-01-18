\name{gen.findDistance}
\alias{gen.findDistance}
\title{Minimum genetic distance between two individuals}
\description{Returns the minimum distance (number of meioses) between the specified individuals through one specified ancestor.}
\usage{gen.findDistance(gen, individuals, ancestor)}
\arguments{
\item{gen}{An object of class GLgen obtained with gen.genealogy, gen.lineages or gen.branching. Required.}
\item{individuals}{A vector of individual id numbers between which to calculate the distance. Required.}
\item{ancestor}{A common ancestor to the specified individuals. Required.}
}
\value{returns a numeric value}
\seealso{

\code{\link{gen.genealogy}}
\code{\link{gen.founder}}
\code{\link{gen.findMRCA}}
\code{\link{gen.findFounders}}
\code{\link{gen.find.Min.Distance.MRCA}}
}
\examples{
data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
gen.findDistance(genJi, individuals=c(1,29), ancestor=17)
gen.findDistance(genJi, individuals=c(1,29), ancestor=26)

data(genea140) 
gen140<-gen.genealogy(genea140) 
gen.findDistance(gen140, individuals=c(409033,408728), ancestor=38714)
gen.findDistance(gen140, individuals=c(408728,408828), ancestor=95080)



}
\keyword{manip}
