\name{gen.implexVar}
\alias{gen.implexVar}
\title{Variance of genealogical implex}
\description{Returns the variance of the implex index across probands.}
\usage{gen.implexVar( gen, pro = "0", onlyNewAnc = "F", genNo = -1, ...)}
\arguments{
\item{gen}{An object of class GLgen obtained with gen.genealogy, gen.lineages or gen.branching. Required.}
\item{pro}{Vector of proband id numbers to be included. Default is 0, which will select all individuals without children.}
\item{onlyNewAnc}{If onlyNewAnc="F" (default), all ancestors will be considered.
	If onlyNewAnc="T", only new ancestors will be counted
	 (i.e., an ancestor is not counted again if it has already been counted in another generation).}
\item{genNo}{Vector of generation numbers at which the implex should be calculated.
	Default is -1, which calculates the implex at each generation.}
\item{...}{Additionnal arguments to be passed to methods.}
}
\value{A data frame with one column and each generation at which the variance is calculated as rows}
\seealso{

\code{\link{gen.gc}}
}
\examples{
data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
gen.implex(genJi, type="IND")
gen.implexVar(genJi)
# For the 5th generation
gen.implexVar(genJi, genNo=5)

\donttest{
data(genea140) 
gen140<-gen.genealogy(genea140) 
gen.implex(gen140)
gen.implex(gen140, pro=c(454422, 676521, 677273, 717537, 717634, 717709, 868572), type="IND")
gen.implexVar(gen140, pro=c(454422, 676521, 677273, 717537, 717634, 717709, 868572), type="IND")
}
}
\references{
Cazes P, Cazes MH. (1996) Comment mesurer la profondeur genealogique d'une ascendance? Population (French Ed) 51:117-140.
}
\keyword{manip}

