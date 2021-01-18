\name{gen.phi}
\alias{gen.phi}
\title{Kinship coefficient}
\description{Returns the kinship coefficients between pairs of individuals}
\usage{gen.phi(gen, pro, depthmin=(gen.depth(gen)-1), depthmax=(gen.depth(gen)-1), MT=F)}
\arguments{
\item{gen}{An object of class GLgen obtained with gen.genealogy, gen.lineages or gen.branching. Required.}
\item{pro}{Vector of proband id numbers to be included. Required.}
\item{depthmin}{Minimum genealogical depth to consider in the calculation. Default is the whole genealogy.}
\item{depthmax}{Maximum genealogical depth to consider in the calculation. Default is the whole genealogy.}
\item{MT}{Allows parallel computing when set to TRUE. Default is MT=F.}
}
\value{A matrix or a GLmultiMatrix object depending on the number of generations treated.
	  GLmultiMatrix is an array of matrices, one for each depth. Array of size ('depthMax'-'depthMin') and 
	  matrices of size 'length(pro)' * 'length(pro)'.
	  The matrix object, also of size 'length(pro)' * 'length(pro)', is returned when ('depthMax'-'depthMin') equals 1.}
\seealso{

\code{\link{gen.genealogy}}
\code{\link{gen.f}}
}
\examples{
data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
kinship<-gen.phi(genJi)
kinship

kinship_allgenerations<-gen.phi(genJi, depthmin =1)
kinship_allgenerations
# 7th generations back in time is equivalent to considering all generations
kinship_allgenerations <- unclass(kinship_allgenerations)
kinship_allgenerations[,,7]==kinship

kinship_allgenerations[1,2,]

# Plot of kinship varying according to number of generations considered
plot(1:7,kinship_allgenerations[1,2,], type="b", xlab="Generation", ylab="Kinship value", 
     ylim=c(0,0.6), pch=0)
points(1:7,kinship_allgenerations[1,3,], type="b", lty=12, pch=1)
legend("topright", legend=c("Individuals 1 and 2", "Individuals 2 and 29"),lty=c(1,12), pch=c(0,1))
}
\references{
Malecot G. (1948) Les mathematiques de l'heredite. Paris: Masson, p 65.

Thompson EA. (1986) Pedigree Analysis in Human Genetics. Baltimore, MD, USA: Johns Hopkins University Press, p 25.

Karigl G. (1981) A recursive algorithm for the calculation of identity coefficients. Ann Hum Genet 45:299-305.
}
\keyword{manip}

