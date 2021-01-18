\name{gen.simuSample}
\alias{gen.simuSample}
\title{Gene dropping simulations - Sample}
\description{Returns the number of alleles transmitted to specified probands from ancestors}
\usage{gen.simuSample(gen, pro, ancestors, stateAncestors, simulNo = 5000)}
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
\value{A matrix with number of columns equal to the number of simulations and number of rows equal to the number of probands.}
\seealso{

\code{\link{gen.genealogy}}
\code{\link{gen.simuProb}}
\code{\link{gen.simuSet}}
\code{\link{gen.simuSampleFreq}}
}
\examples{
data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
simu_1000<-gen.simuSample(genJi, pro=c(1,29), ancestors=c(20,25), stateAncestors=c(2,1),
                          simulNo=1000)
# Number of alleles received by probants
table(simu_1000)
# Number of alleles received by simulation
table(colSums(simu_1000))
# Number of alleles received by each probant
table(simu_1000[1,],simu_1000[2,])

data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
simu_5000<-gen.simuSample(genJi, pro=c(1,29), ancestors=c(20,25), stateAncestors=c(2,1),
                          simulNo=5000)
# Number of alleles received by probants
table(simu_5000)
# Number of alleles received by simulation
table(colSums(simu_5000))
# Number of alleles received by each probant
table(simu_5000[1,],simu_5000[2,])

}
\keyword{datagen}
