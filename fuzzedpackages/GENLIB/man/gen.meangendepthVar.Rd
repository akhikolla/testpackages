\name{gen.meangendepthVar}
\alias{gen.meangendepthVar}
\title{Variance of genealogical depth}
\description{Returns the variance of the genealogical depth}
\usage{gen.meangendepthVar( gen, pro = "0", type = "MEAN", ...)}
\arguments{
\item{gen}{An object of class GLgen obtained with gen.genealogy, gen.lineages or gen.branching. Required.}
\item{pro}{Vector of proband id numbers to be included. Default is 0, which will select all individuals without children.}
\item{type}{If type="MEAN" (default), the average of genealogical depth variances (over all probands) is returned.
	If type="IND", the variance of the genealogical depth is calculated for each specified proband.}
\item{...}{Option to pass additionnal arguments automaticaly between methods. Internal use only.}
}
\value{A data frame with only one numeric value when type is "MEAN".
 	  When type is "IND", the number of rows equals the number of probands specified.}
\seealso{
\code{\link{gen.gc}}
}
\examples{
data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
gen.meangendepth(genJi, type="IND")
gen.meangendepthVar(genJi, type="IND")
gen.meangendepthVar(genJi, type="MEAN")

data(genea140) 
gen140<-gen.genealogy(genea140) 
gen.meangendepth(gen140)
probands <- c(454422,676521,677273,717537,717634,717709,868572)
gen.meangendepth(gen140, pro=probands)
gen.meangendepthVar(gen140, pro=probands)
gen.meangendepthVar(gen140, pro=probands, type="MEAN")

\donttest{
gen.meangendepth(gen140, pro=probands, type="IND")
gen.meangendepthVar(gen140, pro=probands, type="IND")
}
}
\references{
Cazes P, Cazes MH. (1996) Comment mesurer la profondeur genealogique d'une ascendance? Population (French Ed) 51:117-140.

Kouladjian K. (1986) Une mesure d'entropie genealogique. Chicoutimi, SOREP, Document III-C-43.

De Brakaeleer M, Bellis G. (1994) Genealogies et reconstitutions de familles en genetique humaine. Dossiers et Recherches, no 43, INED, Paris.
}
\keyword{manip}

