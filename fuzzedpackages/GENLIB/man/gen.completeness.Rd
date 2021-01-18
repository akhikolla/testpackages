\name{gen.completeness}
\alias{gen.completeness}
\title{Completeness}
\description{Returns the completeness of the genealogical data for the specified probands}
\usage{gen.completeness( gen, pro="0", genNo=-1, type="MEAN", ...)}
\arguments{
\item{gen}{An object of class GLgen obtained with gen.genealogy, gen.lineages or gen.branching. Required.}
\item{pro}{Vector of proband id numbers to be included. Default is 0, which will select all individuals without children.}
\item{genNo}{Vector of generation numbers at which completeness should be calculated.
	Default is -1, which calculates completeness at each generation.}
\item{type}{If type="MEAN" (default), mean completeness over all specified probands is calculated.
	If type="IND", completeness is calculated for each specified proband.}
\item{...}{Option to pass additionnal arguments automaticaly between methods. Internal use only.}
}
\value{A data frame with each generation at which the completeness is calculated as rows and one column when type is "MEAN".
	  When type is "IND", the number of columns equals the number of probands specified.}
\seealso{
\code{\link{gen.genealogy}}
\code{\link{gen.rec}}
\code{\link{gen.occ}}
\code{\link{gen.implex}}
\code{\link{gen.meangendepth}}
}
\examples{
data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
gen.completeness(genJi, type="IND")
# For the 4th generation
gen.completeness(genJi, type="IND", genNo=4)

data(genea140) 
gen140<-gen.genealogy(genea140) 
gen.completeness(gen140)

}
\references{
Cazes P, Cazes MH. (1996) Comment mesurer la profondeur genealogique d'une ascendance? Population (French Ed) 51:117-140.
}
\keyword{manip}

