\name{gen.pro}
\alias{gen.pro}
\title{Get proband id numbers}
\description{Returns the id numbers of the probands. Probands are defined as the individuals without children in the genealogy.}
\usage{gen.pro( gen, ...)}
\arguments{
\item{gen}{An object of class GLgen obtained with gen.genealogy, gen.lineages or gen.branching. Required.}
\item{...}{Option to pass additionnal arguments automaticaly between methods. Internal use only.}
}
\value{returns a vector of integer}
\seealso{
\code{\link{gen.genealogy}}
\code{\link{gen.founder}}
\code{\link{gen.half.founder}}
\code{\link{gen.parent}}
\code{\link{gen.children}}
}
\examples{
data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
# There are 3 probants
gen.pro(genJi)


data(genea140) 
gen140<-gen.genealogy(genea140)
gen.pro(gen140)
# There are 140 probants

}
\keyword{manip}
