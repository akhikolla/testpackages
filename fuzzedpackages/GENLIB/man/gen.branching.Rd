\name{gen.branching}
\alias{gen.branching}
\title{Genealogy subset}
\description{Creates an object of class GLgen from an existing GLgen object by selecting specified individuals.}
\usage{gen.branching( gen, pro = 0, ancestors = gen.founder(gen), bflag = 0)}
\arguments{
\item{gen}{An object of class GLgen obtained with gen.genealogy, gen.lineages or gen.branching. Required.}
\item{pro}{Vector of proband id numbers to be included. Default is 0, which will select all individuals without children.}
\item{ancestors}{Vector of ancestors id numbers to be included. Default includes all founders in the original GLgen object.}
\item{bflag}{If set to 0 (default and recommended), probands and founders are filtered to quicken the calculations.}
}
\value{returns a GLgen object}
\note{Note that if there are no links between some probands and ancestors, these ids will not be included.}
\seealso{
\code{\link{gen.genealogy}}
\code{\link{gen.graph}}
\code{\link{gen.lineages}}
}
\examples{
data(geneaJi) 
genJi<-gen.genealogy(geneaJi)
genJi_part<-gen.branching(genJi,pro=c(2,28)) 

# Plots of original genealogy and of the branched version
layout(matrix(1:2,c(1,2),byrow=TRUE))
gen.graph(genJi) 
mtext("Original", line=2, cex=1.2)
gen.graph(genJi_part) 
mtext("Branched tree\nfor individuals 2 and 28", line=1, cex=1.2)

}
\keyword{manip}
