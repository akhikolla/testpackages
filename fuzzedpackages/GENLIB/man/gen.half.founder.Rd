\name{gen.half.founder}
\alias{gen.half.founder}
\title{Get half-founder id numbers}
\description{Returns the id numbers of the half-founders. Half-founders are defined as the individuals with only one known parent in the genealogy (i.e., either mother id=0 or father id=0).}
\usage{gen.half.founder( gen, ...)}
\arguments{
\item{gen}{An object of class GLgen obtained with gen.genealogy, gen.lineages or gen.branching. Required.}
\item{...}{Option to pass additionnal arguments automaticaly between methods. Internal use only.}
}
\value{returns a vector of integer}
\seealso{
\code{\link{gen.genealogy}}
\code{\link{gen.pro}}
\code{\link{gen.founder}}
\code{\link{gen.parent}}
}
\examples{
data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
# There are 2 half-founders
gen.half.founder(genJi)
}
\keyword{manip}
