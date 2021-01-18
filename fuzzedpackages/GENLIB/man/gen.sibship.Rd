\name{gen.sibship}
\alias{gen.sibship}
\title{Get id numbers of siblings}
\description{Returns the id numbers of the siblings of specified individuals.}
\usage{gen.sibship( gen, individuals, halfSibling = "T", ...)}
\arguments{
\item{gen}{An object of class GLgen obtained with gen.genealogy, gen.lineages or gen.branching. Required.}
\item{individuals}{Vector of individual id numbers. Required.}
\item{halfSibling}{If halfSibling="T" (default) then ids of halfsiblings are also returned. halfSibling="F" returns only full sibling ids.}
\item{...}{Option to pass additionnal arguments automaticaly between methods. Internal use only.}
}
\value{returns a vector of integer}
\seealso{
\code{\link{gen.genealogy}}
\code{\link{gen.pro}}
\code{\link{gen.founder}}
\code{\link{gen.children}}
\code{\link{gen.parent}}
}
\examples{
data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
gen.sibship(genJi, individuals=21,halfSibling=TRUE) 


data(genea140) 
gen140<-gen.genealogy(genea140) 
sibs<-gen.sibship(gen140,individuals=10174, halfSibling=FALSE)
gen.graph(gen140, pro=c(10174,sibs))

}
\keyword{manip}
