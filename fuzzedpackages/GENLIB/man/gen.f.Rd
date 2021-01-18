\name{gen.f}
\alias{gen.f}
\title{Inbreeding coefficient}
\description{Returns the inbreeding coefficients of the specified probands}
\usage{gen.f( gen, pro, depthmin= (gen.depth(gen)-1), depthmax= (gen.depth(gen)-1))}
\arguments{
\item{gen}{An object of class GLgen obtained with gen.genealogy, gen.lineages or gen.branching. Required.}
\item{pro}{Vector of proband id numbers to be included.}
\item{depthmin}{Minimum genealogical depth to consider in the calculation. Default is the whole genealogy.}
\item{depthmax}{Maximum genealogical depth to consider in the calculation. Default is the whole genealogy.}
}
\value{A vector or GLmultiVector object depending on the number of generations treated.
	  GLmultiVector is an array of vectors, one for each depth. Array of size ('depthMax'-'depthMin'), vectors of size 'length(pro)'.
	  A vector of size 'length(pro)' is returned when ('depthMax'-'depthMin') equals 1.}
\seealso{
\code{\link{gen.genealogy}}
\code{\link{gen.phi}}
}
\examples{
data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
f_values<-gen.f(genJi)
f_values
f_allgen<-gen.f(genJi, depthmin=1)  
f_allgen<-unclass(f_allgen)

plot(1:7,f_allgen[1,],type="b",xlab="Generation",ylab="Inbreeding values",ylim=c(0,0.25),pch=0)
points(1:7,f_allgen[3,],  type="b", lty=12, pch=1)
legend("topright", legend=c("Individual 1", "Individual 29"),lty=c(1,12), pch=c(0,1))

\donttest{
data(genea140) 
gen140<-gen.genealogy(genea140)
inbreeding_val<-gen.f(gen140)
boxplot(inbreeding_val, horizontal=TRUE, xlab="Inbreeding values")
}
}
\references{
Malecot G. (1948) Les mathematiques de l'heredite. Paris: Masson, p 65.

Thompson EA. (1986) Pedigree Analysis in Human Genetics. Baltimore, MD, USA: Johns Hopkins University Press, p 25.

Karigl G. (1981) A recursive algorithm for the calculation of identity coefficients. Ann Hum Genet 45:299-305.
}
\keyword{manip}

