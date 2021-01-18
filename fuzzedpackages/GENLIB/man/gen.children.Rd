\name{gen.children}
\alias{gen.children}
\title{Get id numbers of children}
\description{Returns the id numbers of the children of specified individuals.}
\usage{gen.children( gen, individuals, ...)}
\arguments{
\item{gen}{An object of class GLgen obtained with gen.genealogy, gen.lineages or gen.branching. Required.}
\item{individuals}{Vector of individual id numbers. Required.}
\item{...}{Option to pass additionnal arguments automaticaly between methods. Internal use only.}
}
\value{returns a vector of integer}
\seealso{
\code{\link{gen.genealogy}}
\code{\link{gen.pro}}
\code{\link{gen.founder}}
\code{\link{gen.sibship}}
\code{\link{gen.parent}}
}
\examples{
data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
gen.children(genJi, individuals=14) 

data(genea140) 
gen140<-gen.genealogy(genea140) 
children51052<-gen.children(gen140,51052)
gen.graph(gen140, pro=c(children51052))
}
\keyword{manip}
