\name{gen.findFounders}
\alias{gen.findFounders}
\title{Find common founder ancestors}
\description{Returns all the ancestors that are founders for specified individuals in the genealogy.}
\usage{gen.findFounders(gen, individuals, NbProcess=parallel::detectCores()-1)}
\arguments{
\item{gen}{An object of class GLgen obtained with gen.genealogy, gen.lineages or gen.branching. Required.}
\item{individuals}{Vector of individual id numbers for which to find the founders. Required}
\item{NbProcess}{Number of processes to use when running this function. Default=parallel::detectCores()-1}
}
\value{A vector of integers}
\note{Uses slave processes to make the search faster. Those slave processes launch the function gen.getFoundersMPI.}
\seealso{
\code{\link{gen.genealogy}}
\code{\link{gen.founder}}
\code{\link{gen.findMRCA}}
\code{\link{gen.findDistance}}
\code{\link{gen.find.Min.Distance.MRCA}}
}
\examples{
data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
gen.founder(genJi)
gen.findFounders(genJi, individuals=c(1,29), NbProcess = 1)

data(genea140) 
gen140<-gen.genealogy(genea140) 
gen2<-gen.branching(gen140,pro=c(409033,408728)) 
# 615 founders in genealogy of #409033 and #408728
all_founders<-gen.founder(gen2) 
length(all_founders)

\donttest{
all_commonFounders<-gen.findFounders(gen140, individuals=c(409033,408728), NbProcess = 1)
length(all_commonFounders) # 127 founders common to #409033 and #408728
}
}
\keyword{manip}
