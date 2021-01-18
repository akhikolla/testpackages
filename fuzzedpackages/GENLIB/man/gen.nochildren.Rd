\name{gen.nochildren}
\alias{gen.nochildren}
\title{Number of children}
\description{Returns the number of children for specified individuals}
\usage{gen.nochildren( gen, individuals)}
\arguments{
\item{gen}{An object of class GLgen obtained with gen.genealogy, gen.lineages or gen.branching. Required.}
\item{individuals}{A vector of ids of the individuals on each of whom the number of chidren is to be returned. Required.}
}
\value{A vector containing the number of children for each individual sepcified.}
\seealso{
\code{\link{gen.genealogy}}
\code{\link{gen.nowomen}}
\code{\link{gen.noind}}
\code{\link{gen.nomen}}
\code{\link{gen.children}}
}
\examples{
data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
gen.children(genJi,14)
gen.nochildren(genJi,individuals=c(1,12,14,20))

data(genea140) 
gen140<-gen.genealogy(genea140) 
children51052<-gen.children(gen140,51052)
gen.nochildren(gen140,51052)
gen.graph(gen140, pro=children51052)
}
\keyword{manip}
