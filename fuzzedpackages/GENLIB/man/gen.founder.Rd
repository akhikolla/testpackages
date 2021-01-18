\name{gen.founder}
\alias{gen.founder}
\title{Get founder id numbers}
\description{Returns the id numbers of the founders. Founders are defined as the individuals without parents in the genealogy (i.e., mother id=0, father id=0).}
\usage{gen.founder( gen, ...)}
\arguments{
\item{gen}{An object of class GLgen obtained with gen.genealogy, gen.lineages or gen.branching. Required.}
\item{...}{Option to pass additionnal arguments automaticaly between methods. Internal use only.}
}
\value{A vector of integer}
\seealso{
\code{\link{gen.genealogy}}
\code{\link{gen.pro}}
\code{\link{gen.half.founder}}
\code{\link{gen.parent}}
}
\examples{
data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
# There are 6 founders
gen.founder(genJi)

data(genea140) 
gen140<-gen.genealogy(genea140)
founder140<-gen.founder(gen140)  
# There are 7399 founders
length(founder140)

}
\keyword{manip}
