\name{gen.simuSampleFreq}
\alias{gen.simuSampleFreq}
\title{Gene dropping simulations - Frequencies}
\description{Returns the number of alleles transmitted to specified probands from ancestors}
\usage{gen.simuSampleFreq( gen, pro, ancestors, stateAncestors, simulNo = "5000")}
\arguments{
\item{gen}{An object of class GLgen obtained with gen.genealogy, gen.lineages or gen.branching. Required.}
\item{pro}{Vector of proband id numbers to be included. Required.}
\item{ancestors}{Vector of ancestors id numbers to be included. Required.}
\item{stateAncestors}{Required vector indicating, for each ancestor in ancestors, the genotype state:

	0 = no disease allele present

	1 = 1 disease allele present

	2 = 2 disease alleles present}
\item{simulNo}{Number of simulations to perform. Default is 5000.}
}
\value{A data frame containing for each row (one row per proband):
\strong{Alleles.transmitted.0} 
\strong{Alleles.transmitted.1} 
\strong{Alleles.transmitted.2} 
}
\seealso{

\code{\link{gen.genealogy}}
\code{\link{gen.simuProb}}
\code{\link{gen.simuSet}}
\code{\link{gen.simuSample}}
}
\examples{
data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
gen.simuSampleFreq(genJi, pro=c(1,29), ancestors=c(20,25),stateAncestors=c(2,1), simulNo = 1000)

}
\keyword{datagen}
