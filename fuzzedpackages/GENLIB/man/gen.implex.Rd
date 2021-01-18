\name{gen.implex}
\alias{gen.implex}
\title{Genealogical implex}
\description{Returns the genealogical implex index (a measure of pedigree collapsing) for the specified probands.}
\usage{gen.implex( gen, pro = "0", genNo = -1, type = "MEAN", onlyNewAnc = "F", ...)}
\arguments{
\item{gen}{An object of class GLgen obtained with gen.genealogy, gen.lineages or gen.branching. Required.}
\item{pro}{Vector of proband id numbers to be included. Default is 0, which will select all individuals without children.}
\item{genNo}{Vector of generation numbers at which the implex should be calculated.
	Default is -1, which calculates the implex at each generation.}
\item{type}{If type="MEAN" (default), implex index values are averaged over all specified probands.
	If type="IND", the implex index is calculated for all specified proband together.
	If type="ALL", the implex index is calculated for each specified proband.}
\item{onlyNewAnc}{If onlyNewAnc="F" (default), all ancestors will be considered.
	If onlyNewAnc="T", only new ancestors will be counted
	 (i.e., an ancestor is not counted again if it has already been counted in another generation).}
\item{...}{Option to pass additionnal arguments automaticaly between methods. Internal use only.}
}
\value{A data frame with each generation at which the completeness is calculated as rows and one column when type is "MEAN".
 	  When type is "IND", the number of columns equals the number of probands specified.}
\seealso{

\code{\link{gen.genealogy}}
\code{\link{gen.occ}}
\code{\link{gen.rec}}
\code{\link{gen.meangendepth}}
\code{\link{gen.completeness}}
\code{\link{gen.gc}}
}
\examples{
data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
gen.implex(genJi)
gen.implex(genJi, type="IND")
# For the 5th generation
gen.implex(genJi, type="IND", genNo=5) 

\donttest{
data(genea140) 
gen140<-gen.genealogy(genea140) 
gen.implex(gen140)
gen.implex(gen140, pro=c(454422, 676521, 677273, 717537, 717634, 717709, 868572))
gen.implex(gen140, pro=c(454422, 676521, 677273, 717537, 717634, 717709, 868572), type="IND")
}
}
\references{
Cazes P, Cazes MH. (1996) Comment mesurer la profondeur genealogique d'une ascendance? Population (French Ed) 51:117-140.
}
\keyword{manip}

