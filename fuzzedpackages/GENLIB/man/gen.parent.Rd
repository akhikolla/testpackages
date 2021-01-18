\name{gen.parent}
\alias{gen.parent}
\title{Get id numbers of parents}
\description{Returns the id numbers of the parents of specified individuals.}
\usage{gen.parent( gen, individuals, output = "FaMo", ...)}
\arguments{
\item{gen}{An object of class GLgen obtained with gen.genealogy, gen.lineages or gen.branching. Required.}
\item{individuals}{Vector of individual id numbers. Required.}
\item{output}{If output="FaMo" (default) then both mothers and fathers are included. "Mo" outputs mothers only and "Fa", fathers only.}
\item{...}{Option to pass additionnal arguments automaticaly between methods. Internal use only.}
}
\value{returns a list containing the following:
\strong{Fathers} 
\strong{Mothers} 
}
\seealso{
\code{\link{gen.genealogy}}
\code{\link{gen.pro}}
\code{\link{gen.founder}}
\code{\link{gen.children}}
\code{\link{gen.sibship}}
}
\examples{
data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
gen.parent(genJi, individuals=c(3,21,29)) 

}
\keyword{manip}
