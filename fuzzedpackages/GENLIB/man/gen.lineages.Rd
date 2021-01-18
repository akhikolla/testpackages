\name{gen.lineages}
\alias{gen.lineages}
\title{Create object of class GLgen for maternal or paternal lineages}
\description{Creates an object of class GLgen that contains maternal or paternal lineages selected from input data in pedigree format}
\usage{gen.lineages(ped, pro = "0", maternal = "T", ...)}
\arguments{
\item{ped}{A matrix or data frame with the following columns: individual id, father id, mother id, sex. Required.
	Individal id numbers must be numeric and unique.
	If an individual does not have a father and/or mother, the father and/or mother id numbers must be set to 0.
	All non-zero father and mother id numbers must also appear in the individual id column.
	Sex must be equal to M or 1 for males and F or 2 for females.
	The sex column is optional for this function but necessary for some other functions using GLgen objects.}
\item{pro}{Vector of individual id numbers for which lineages should be included. Optional.}
\item{maternal}{If mat="T" (default), maternal lineages are selected. mat="F" returns paternal lineages.}
\item{...}{Option to pass additionnal arguments automaticaly between methods. Internal use only.}
}
\value{returns a GLgen object}
\seealso{
\code{\link{gen.genealogy}}
\code{\link{gen.graph}}
\code{\link{gen.branching}}
\code{\link{gen.genout}}
}
\examples{
data(geneaJi) 
genJi <- gen.genealogy(geneaJi)
genJi_MaLi<-gen.lineages(geneaJi, maternal = TRUE) 
genJi_FaLi<-gen.lineages(geneaJi, maternal = FALSE) 

# Plots of original genealogy and maternal and paternal lineages
layout(matrix(1:3,c(1,3),byrow=TRUE), widths =c(3,1,1), heights = 1)
gen.graph(genJi) 
mtext("Original", line=2)
gen.graph(genJi_MaLi) 
mtext("Maternal\nlineages", line=1)
gen.graph(genJi_FaLi) 
mtext("Paternal\nlineages", line=1)

}
\keyword{manip}
