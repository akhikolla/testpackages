\name{gen.mean}
\alias{gen.mean}
\title{Mean number of generations}
\description{Returns the average number of generations between all probands and the individuals specified. Probands are defined as the individuals without children in the genealogy.}
\usage{gen.mean( gen, individuals)}
\arguments{
\item{gen}{An object of class GLgen obtained with gen.genealogy, gen.lineages or gen.branching. Required.}
\item{individuals}{A vector of ids specifying the individuals to include in the calculation. Required.}
}
\value{returns a numeric value}
\seealso{

\code{\link{gen.genealogy}}
\code{\link{gen.min}}
\code{\link{gen.max}}
\code{\link{gen.depth}}
\code{\link{gen.meangendepth}}
}
\examples{
data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
gen.min(genJi,c(17,26))
gen.mean(genJi,c(17,26))
gen.max(genJi,c(17,26))


data(genea140) 
gen140<-gen.genealogy(genea140) 
gen.min(gen140,c(18311,18430))
gen.mean(gen140,c(18311,18430))
gen.max(gen140,c(18311,18430))
}
\keyword{manip}
