\name{gen.occ}
\alias{gen.occ}
\title{Ancestor occurences}
\description{Returns the number of times that the specified ancestors are present in the genalogies of the specified probands.}
\usage{gen.occ( gen, pro = "0", ancestors = "0", typeOcc = "IND", ...)}
\arguments{
\item{gen}{An object of class GLgen obtained with gen.genealogy, gen.lineages or gen.branching. Required.}
\item{pro}{Vector of proband id numbers to be included. Default is 0, which will select all individuals without children.}
\item{ancestors}{Vector of ancestors id numbers to be included. Default is 0, which will select all individuals without parents.}
\item{typeOcc}{If typeOcc="IND" (default), the number of occurences for each proband will be returned.
	If typeOcc="TOTAL", the sum of all occurences over all probands will be returned.}
\item{...}{Option to pass additionnal arguments automaticaly between methods. Internal use only.}
}
\value{A matrix with number of lines equal to the number of ancestors included and the number of columns equal to the number of probands included if typeOcc="BRUT" or only one column if typeOcc="TOTAL"}
\seealso{

\code{\link{gen.genealogy}}
\code{\link{gen.rec}}
\code{\link{gen.implex}}
\code{\link{gen.meangendepth}}
\code{\link{gen.gc}}
}
\examples{
data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
# Number of occurences of ancestors in the specified proband's genealogy
gen.occ(genJi, pro=c(1,29), ancestors=c(3,6,10,12,14,16,18,20,26,28))

data(genea140) 
gen140<-gen.genealogy(genea140) 
gen.occ(gen140, pro=c(454422,676521,677273,717537,717634,717709,868572), 
        ancestors=c(18305,18528,31114,18325))

}
\keyword{manip}
