\name{gen.depth}
\alias{gen.depth}
\title{Genealogy depth}
\description{Returns the number of generations in the genealogy from a GLgen object.}
\usage{gen.depth(gen)}
\arguments{
\item{gen}{An object of class GLgen obtained with gen.genealogy, gen.lineages or gen.branching. Required.}
}
\value{returns a vector of integer}
\seealso{
\code{\link{gen.genealogy}}
\code{\link{gen.min}}
\code{\link{gen.mean}}
\code{\link{gen.min}}
\code{\link{gen.meangendepth}}
}
\examples{
data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
gen.depth(genJi)


data(genea140) 
gen140<-gen.genealogy(genea140) 
gen.depth(gen140)

}
\keyword{attribute}
