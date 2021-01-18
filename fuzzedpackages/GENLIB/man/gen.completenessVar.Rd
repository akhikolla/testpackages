\name{gen.completenessVar}
\alias{gen.completenessVar}
\title{Variance of completeness index}
\description{Returns the variance of the completeness index across probands.}
\usage{gen.completenessVar( gen, pro = "0", genNo = -1, ...)}
\arguments{
\item{gen}{An object of class GLgen obtained with gen.genealogy, gen.lineages or gen.branching. Required.}
\item{pro}{Vector of proband id numbers to be included. Default is 0, which will select all individuals without children.}
\item{genNo}{Vector of generation numbers at which completeness should be calculated.
	Default is -1, which calculates completeness at each generation.}
\item{...}{Option to pass additionnal arguments automaticaly between methods. Internal use only.}
}
\value{A data frame with one column and each generation at which the variance is calculated as rows.}
\seealso{
\code{\link{gen.gc}}
}
\examples{
data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
gen.completeness(genJi, type="IND")
gen.completenessVar(genJi)
# For the 4th generation
gen.completenessVar(genJi, genNo=4)

\donttest{
data(genea140) 
gen140<-gen.genealogy(genea140) 
gen.completeness(gen140)
gen.completenessVar(gen140)
}
}
\references{
Cazes P, Cazes MH. (1996) Comment mesurer la profondeur genealogique d'une ascendance? Population (French Ed) 51:117-140.
}
\keyword{manip}

